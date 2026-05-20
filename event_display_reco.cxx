// =============================================================================
//  event_display_reco.cxx
//
//  REve-based event display showing the pi0 reconstruction:
//
//    * Loads a GDML geometry (same filtered-volumes approach as the original
//      event_display.cxx).
//    * Reads three data sources, all aligned by event index:
//        - "calo_events" (truth-hit cloud, same as the original);
//        - "truth_particles" (truth pi0 vertex + photon directions);
//        - "cluster_tree" (reconstructed clusters and matches).
//    * Picks the event with the smallest |z_reco_POCA - z_truth_vertex|.
//    * For the chosen event, displays:
//        - the original edep point cloud as backdrop (5 log-binned colors);
//        - two cylinders, one per reconstructed photon, with axis along the
//          reco line, length matched to the cluster hits' z extent, and
//          radius equal to the transverse MIP-weighted RMS of the cluster
//          hits about the fitted line;
//        - two SOLID lines pointing from each cluster cylinder centre to
//          the reconstructed (POCA) vertex;
//        - two DASHED lines pointing from each cluster cylinder centre to
//          the truth vertex;
//        - two marker points: green star = truth vertex, magenta cross =
//          reco POCA vertex.
//
//  Build with the same CMakeLists.txt: add this file as a new executable
//  target alongside event_display.cxx.
//
//  Run:   ./event_display_reco <detector.gdml> <truth+calo.root> <cluster.root>
//                              [--event N] [--tree-truth calo_events]
//                              [--tree-particles truth_particles]
//                              [--tree-cluster cluster_tree]
// =============================================================================

#include <ROOT/REveManager.hxx>
#include <ROOT/REveScene.hxx>
#include <ROOT/REveViewer.hxx>
#include <ROOT/REveElement.hxx>
#include <ROOT/REveGeoShape.hxx>
#include <ROOT/REveBoxSet.hxx>
#include <ROOT/REvePointSet.hxx>
#include <ROOT/REveStraightLineSet.hxx>
#include <ROOT/REveRGBAPalette.hxx>
#include <ROOT/REveTrans.hxx>

#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TGeoVolume.h>
#include <TGeoMatrix.h>
#include <TGeoShape.h>
#include <TGeoBBox.h>
#include <TGeoTube.h>
#include <TGeoSphere.h>

#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TColor.h>
#include <TString.h>
#include <Math/Vector3D.h>

#include <iostream>
#include <limits>
#include <vector>
#include <regex>
#include <string>
#include <algorithm>
#include <cmath>
#include <map>

namespace REX = ROOT::Experimental;

// =============================================================================
//  Reconstruction helpers (mirroring cluster_pca.C / reco_compare.py)
// =============================================================================
struct Vec3 { double x, y, z; };
struct Line { Vec3 p, d; bool ok = false; };  // d is a unit direction vector

struct LineFit2D {
    double slope = 0.0;
    double intercept = 0.0;
    double zRef = 0.0;
    bool   ok = false;
};

// Weighted line fit p = slope*z + (intercept at z = zRef)
static LineFit2D fitLine2D(const std::vector<double>& z,
                           const std::vector<double>& p,
                           const std::vector<double>& w)
{
    LineFit2D F;
    if (z.size() < 3) return F;
    double sw = 0, swz = 0, swp = 0, swzz = 0, swzp = 0;
    for (size_t i = 0; i < z.size(); ++i) {
        sw   += w[i];
        swz  += w[i] * z[i];
        swp  += w[i] * p[i];
        swzz += w[i] * z[i] * z[i];
        swzp += w[i] * z[i] * p[i];
    }
    if (sw <= 0) return F;
    F.zRef = swz / sw;
    double Szz = swzz - swz * swz / sw;
    double Szp = swzp - swz * swp / sw;
    if (std::abs(Szz) < 1e-9) return F;
    F.slope     = Szp / Szz;
    F.intercept = swp / sw;   // value at z = zRef
    F.ok        = true;
    return F;
}

static LineFit2D iterateHaloCut2D(const std::vector<double>& Z,
                                  const std::vector<double>& P,
                                  const std::vector<double>& W,
                                  double rMaxMm, int maxIter, double convTolMm)
{
    LineFit2D F = fitLine2D(Z, P, W);
    if (!F.ok) return F;
    for (int it = 0; it < maxIter; ++it) {
        std::vector<double> Zc, Pc, Wc;
        for (size_t i = 0; i < Z.size(); ++i) {
            double pred = F.intercept + F.slope * (Z[i] - F.zRef);
            if (std::abs(P[i] - pred) < rMaxMm) {
                Zc.push_back(Z[i]); Pc.push_back(P[i]); Wc.push_back(W[i]);
            }
        }
        if (Zc.size() < 3) break;
        LineFit2D Fnew = fitLine2D(Zc, Pc, Wc);
        if (!Fnew.ok) break;
        bool conv = std::abs(Fnew.slope - F.slope) < convTolMm / 1e4
                 && std::abs(Fnew.intercept - F.intercept) < convTolMm;
        F = Fnew;
        if (conv) break;
    }
    return F;
}

// Combine two 2D fits (z,x) and (z,y) into a 3D line.
static Line combineFits(const LineFit2D& Fx, const LineFit2D& Fy)
{
    Line L;
    if (!Fx.ok || !Fy.ok) return L;
    const double zp = 0.5 * (Fx.zRef + Fy.zRef);
    L.p.x = Fx.intercept + Fx.slope * (zp - Fx.zRef);
    L.p.y = Fy.intercept + Fy.slope * (zp - Fy.zRef);
    L.p.z = zp;
    double dx = Fx.slope, dy = Fy.slope, dz = 1.0;
    double n  = std::sqrt(dx*dx + dy*dy + dz*dz);
    L.d.x = dx / n; L.d.y = dy / n; L.d.z = dz / n;
    L.ok = true;
    return L;
}

// POCA between two 3D lines; returns midpoint and ok flag.
struct PocaResult { Vec3 vtx; double dist; bool ok = false; };
static PocaResult poca(const Line& A, const Line& B)
{
    PocaResult R;
    if (!A.ok || !B.ok) return R;
    Vec3 dr{ B.p.x - A.p.x, B.p.y - A.p.y, B.p.z - A.p.z };
    double b = A.d.x*B.d.x + A.d.y*B.d.y + A.d.z*B.d.z;
    double dnum = A.d.x*dr.x + A.d.y*dr.y + A.d.z*dr.z;
    double enum_ = B.d.x*dr.x + B.d.y*dr.y + B.d.z*dr.z;
    double det = 1.0 - b * b;
    if (std::abs(det) < 1e-12) return R;
    double t1 = (dnum - b * enum_) / det;
    double t2 = (b * dnum - enum_) / det;
    Vec3 P1{ A.p.x + t1*A.d.x, A.p.y + t1*A.d.y, A.p.z + t1*A.d.z };
    Vec3 P2{ B.p.x + t2*B.d.x, B.p.y + t2*B.d.y, B.p.z + t2*B.d.z };
    R.vtx = { 0.5*(P1.x+P2.x), 0.5*(P1.y+P2.y), 0.5*(P1.z+P2.z) };
    double dx = P1.x-P2.x, dy = P1.y-P2.y, dz = P1.z-P2.z;
    R.dist = std::sqrt(dx*dx + dy*dy + dz*dz);
    R.ok = true;
    return R;
}

// One reconstructed photon's cluster info
struct PhotonReco {
    Line   line;
    Vec3   centroid;        // MIP-weighted (x, y, z) of cluster hits
    double z_extent_half;   // half-length along line (~ shower depth/2)
    double r_rms;           // MIP-weighted transverse RMS about the line
    int    hClId = -1;      // H cluster ID
    int    vClId = -1;      // V cluster ID
    bool   ok = false;
};

// =============================================================================
//  EventDisplay (extends the original to consume truth + cluster trees)
// =============================================================================
class EventDisplay {
public:
    EventDisplay()  = default;
    ~EventDisplay() = default;

    void Init();
    void LoadGeometry(const std::string& gdml_path,
                      const std::vector<std::string>& include_regex,
                      int max_depth = -1);

    void OpenEventFile(const std::string& root_path,
                       const std::string& calo_tree    = "calo_events",
                       const std::string& truth_tree   = "truth_particles");
    void OpenClusterFile(const std::string& cluster_path,
                         const std::string& tree_name  = "cluster_tree");

    Long64_t NumEvents() const {
        if (fTreeClus) return fTreeClus->GetEntries();
        if (fTreeCalo) return fTreeCalo->GetEntries();
        if (fTreeTrue) return fTreeTrue->GetEntries();
        return 0;
    }

    // Pick the event with the smallest |z_reco - z_truth| over all events.
    Long64_t PickBestEvent();

    void GotoEvent(Long64_t i);
    REX::REveManager* Manager() { return fEve; }

    void SetHitScale(double s) { fHitScale = s; }
    void SetSceneOffset(double x, double y, double z) {
        fOff[0] = x; fOff[1] = y; fOff[2] = z; fOffValid = true;
    }
    void AutoCenterOnFirstEvent();
    void AutoCenterGeometry();

    // Create a second viewer bound ONLY to the geometry (global) scene, so
    // its camera auto-frames the detector bounding box -- i.e. the camera
    // is "centered and zoomed on the geometry" by default.  Call after the
    // geometry has been loaded and placed.
    void SetupViewers();

    // Diagnostic: print the scene-space (cm) bounding boxes of the
    // geometry, the calo edep hits, the cluster hits, and the truth/reco
    // vertices for event `i`, so frame-alignment problems are visible.
    void DumpScenePlacement(Long64_t i);

private:
    // ---- Geometry loading -------------------------------------------------
    void BuildGeoShapes(TGeoNode* node, const TGeoHMatrix& parent_mtx,
                        const std::vector<std::regex>& patterns,
                        int depth, int max_depth,
                        REX::REveElement* parent);

    // ---- Per-event reconstruction -----------------------------------------
    // Run the cluster_pca algorithm on the *current* cluster_tree entry
    // (which must already be loaded via fTreeClus->GetEntry()).
    // Returns true on success; fills phs[2] and pocaVtx.
    bool ReconstructCurrentEvent(PhotonReco phs[2], Vec3& pocaVtx) const;
    bool TruthFromCurrentEvent(Vec3& truthVtx, Vec3 truthDir[2]) const;

    // ---- REve elements ----------------------------------------------------
    REX::REveManager* fEve         = nullptr;
    REX::REveScene*   fGeoScene    = nullptr;
    REX::REveScene*   fEventScene  = nullptr;
    REX::REveElement* fGeoHolder   = nullptr;
    REX::REveElement* fEventHolder = nullptr;

    // ---- File / tree pointers ---------------------------------------------
    TFile* fFileCalo = nullptr;
    TFile* fFileClus = nullptr;
    TTree* fTreeCalo = nullptr;   // calo_events
    TTree* fTreeTrue = nullptr;   // truth_particles
    TTree* fTreeClus = nullptr;   // cluster_tree

    // Uniform scene scale: scene units == analysis-frame mm * fHitScale.
    // fHitScale is a SIMILARITY transform -- every distance, size and angle
    // is preserved exactly; only the absolute numbers shrink.  It is set
    // below 1 because ROOT 6.38's web REve gives no camera clip-plane
    // control, and a perspective camera cannot resolve the true ~58 m
    // scene span (the depth buffer loses precision -> blank / banded 3D
    // view).  At 0.01 the span becomes ~0.58 m and the perspective 3D view
    // renders correctly.  Applied to hits/cylinders/vertices via toScene()
    // and to the geometry in AutoCenterGeometry() -- all with the SAME
    // factor, so nothing is distorted.
    double fHitScale = 0.01;      // mm -> scaled scene units

    // Rigid scene translation, in SCALED units: subtracted from every
    // world coordinate so the detector bounding-box centre lands on the
    // scene origin (0,0,0).  A pure change of origin -- no distortion.
    // Computed in AutoCenterGeometry from the detector bbox; applied
    // identically to hits, cluster hits, cylinders, vertices (via toScene)
    // and to the geometry itself.  DumpScenePlacement divides fHitScale
    // and fOff back out, so the diagnostic still reports true mm.
    double fOff[3]   = {0, 0, 0};
    bool   fOffValid = false;     // set true once the detector bbox centre
                                  // has been computed
    double fGeoOff[3] = {0, 0, 0};
    bool   fGeoOffValid = false;

    // ---- Branch buffers: calo_events --------------------------------------
    std::vector<double>* fEdep    = nullptr;
    std::vector<double>* fXg      = nullptr;
    std::vector<double>* fYg      = nullptr;
    std::vector<double>* fZg      = nullptr;

    // ---- Branch buffers: truth_particles ---------------------------------
    std::vector<double>* fVx_mm    = nullptr;
    std::vector<double>* fVy_mm    = nullptr;
    std::vector<double>* fVz_mm    = nullptr;
    std::vector<double>* fPx_MeV   = nullptr;
    std::vector<double>* fPy_MeV   = nullptr;
    std::vector<double>* fPz_MeV   = nullptr;
    std::vector<int>*    fFirstGen = nullptr;

    // ---- Branch buffers: cluster_tree -------------------------------------
    std::vector<int>*    fHitOri    = nullptr;
    std::vector<double>* fHitPos    = nullptr;
    std::vector<double>* fHitMip    = nullptr;
    std::vector<double>* fHitZ      = nullptr;
    std::vector<double>* fHitWeight = nullptr;
    std::vector<int>*    fHitClust  = nullptr;
    std::vector<int>*    fClId      = nullptr;
    std::vector<int>*    fClOri     = nullptr;
    std::vector<double>* fClMip     = nullptr;
    std::vector<int>*    fMatchHcl  = nullptr;
    std::vector<int>*    fMatchVcl  = nullptr;
    std::vector<double>* fMatchScore = nullptr;
};

// =============================================================================
//  Boilerplate (mostly cribbed from event_display.cxx)
// =============================================================================
void EventDisplay::Init()
{
    fEve         = REX::REveManager::Create();
    fGeoScene    = fEve->GetGlobalScene();
    fEventScene  = fEve->GetEventScene();
    fGeoHolder   = new REX::REveElement("Detector");
    fGeoScene->AddElement(fGeoHolder);
    fEventHolder = new REX::REveElement("Event");
    fEventScene->AddElement(fEventHolder);
}

// Configure viewers for a FIXED-TARGET calorimeter display.
//
// Why not REveCalo3D/REveCalo2D: those classes model a COLLIDER
// calorimeter -- REveCaloData stores every cell as (etaMin,etaMax,
// phiMin,phiMax) and draws projective towers around a beam axis.  In a
// fixed-target setup the calorimeter is a flat Cartesian cell grid with
// no natural eta/phi, so those classes do not apply.  The hits/cells are
// therefore drawn at their real Cartesian (x,y,z) -- the correct
// representation here.
//
// ROOT 6.38's web REve exposes no camera position/zoom/clip setter
// (REveCamera is "internal"); the only camera control is
// REveViewer::SetCameraType().  The enum (from REveViewer.hxx) includes
// orthographic modes:
//     kCameraOrthoXOY  -- looking down  Z : X horizontal, Y vertical
//     kCameraOrthoXOZ  -- looking along Y : X horizontal, Z vertical
// These are the natural "2D" views for a fixed-target detector:
//   * XOY = beam's-eye view of the calorimeter face,
//   * XOZ = side view, showing depth and the photon directions.
// Orthographic cameras also have no perspective foreshortening and a flat
// depth treatment, which avoids the depth-buffer banding seen when the
// perspective camera is zoomed out over the large (~58 m) scene span.
//
// Three viewers are created, each showing BOTH the geometry and the event
// scene (detector + hits together), differing only in camera:
//   "3D view"   -- perspective       (kCameraPerspXOZ)   : top-left, large
//   "Side XZ"   -- orthographic Y    (kCameraOrthoXOZ)   : top-right
//   "Beam YZ"   -- orthographic X    (kCameraOrthoZOY)   : long bottom strip
// The pre-existing default viewer is left as an extra perspective view.
// (Window placement is done by the EVE client -- drag/resize panels there;
// the layout can be saved client-side.)
void EventDisplay::SetupViewers()
{
    if (!fEve) return;

    // Helper: spawn a viewer, attach both scenes, set its camera type.
    auto makeViewer = [&](const char* name, const char* title,
                          REX::REveViewer::ECameraType cam) {
        REX::REveViewer* v = fEve->SpawnNewViewer(name, title);
        if (!v) {
            std::cerr << "[ED] SetupViewers: could not spawn '"
                      << name << "'\n";
            return;
        }
        // Show detector AND event overlays together in every viewer.
        v->AddScene(fGeoScene);
        v->AddScene(fEventScene);
        v->SetCameraType(cam);
        std::cout << "[ED] spawned viewer '" << name << "'\n";
    };

    // 3D perspective view (top-left, large).
    makeViewer("3D view", "Perspective 3D view",
               REX::REveViewer::kCameraPerspXOZ);
    // Orthographic side elevation (top-right): looking along Y,
    // X horizontal, Z vertical.
    makeViewer("Side XZ", "Orthographic side elevation",
               REX::REveViewer::kCameraOrthoXOZ);
    // Orthographic YZ view (long bottom strip): looking along X, with the
    // beam axis Z horizontal and Y vertical.  This is the view that shows
    // the vertex-to-detector standoff along the beam -- a wide, short
    // viewer suits it.  Orthographic, so no perspective depth banding.
    makeViewer("Beam YZ", "Orthographic YZ view along the beam axis",
               REX::REveViewer::kCameraOrthoZOY);

    // Leave the pre-existing default viewer as an extra perspective view.
    if (auto* def = fEve->GetDefaultViewer()) {
        def->SetCameraType(REX::REveViewer::kCameraPerspXOZ);
    }

    // Reset all viewer cameras to fit their (final) scene contents.
    // RepaintAllViewers lives on REveViewerList (reached via GetViewers()).
    if (auto* vl = fEve->GetViewers()) {
        vl->RepaintAllViewers(/*resetCameras*/ true, /*dropLogicals*/ false);
    }
    std::cout << "[ED] viewers configured (3D perspective + "
                 "orthographic XY beam's-eye + orthographic XZ side)\n";
}

void EventDisplay::LoadGeometry(const std::string& gdml_path,
                                const std::vector<std::string>& include_regex,
                                int max_depth)
{
    TGeoManager* geo = TGeoManager::Import(gdml_path.c_str());
    if (!geo) {
        std::cerr << "[ED] failed to import GDML: " << gdml_path << "\n";
        return;
    }
    std::vector<std::regex> patterns;
    for (const auto& s : include_regex)
        patterns.emplace_back(s, std::regex::ECMAScript);

    TGeoHMatrix id;
    BuildGeoShapes(geo->GetTopNode(), id, patterns, 0, max_depth, fGeoHolder);
    std::cout << "[ED] geometry: " << fGeoHolder->NumChildren() << " shapes\n";
}

void EventDisplay::BuildGeoShapes(TGeoNode* node, const TGeoHMatrix& parent_mtx,
                                  const std::vector<std::regex>& patterns,
                                  int depth, int max_depth,
                                  REX::REveElement* parent)
{
    if (max_depth >= 0 && depth > max_depth) return;
    if (!node) return;
    TGeoHMatrix global = parent_mtx;
    global.Multiply(node->GetMatrix());

    const std::string vname = node->GetVolume()->GetName();
    bool matched = false;
    for (const auto& p : patterns) {
        if (std::regex_search(vname, p)) { matched = true; break; }
    }
    if (matched) {
        TGeoShape* shape = node->GetVolume()->GetShape();
        if (shape) {
            auto* eshape = new REX::REveGeoShape(vname.c_str());
            eshape->SetShape(shape);
            eshape->RefMainTrans().SetFrom(global);
            // Color by volume name pattern.  The full list of materials:
            //   Lead absorber  -> dark gray (the bulk of the calo mass)
            //   WidePVT        -> orange (the wide scintillator bars)
            //   ThinPS         -> yellow (the thin scintillator bars)
            //   HCAL_*         -> green (HCAL scintillator)
            //   IronPlate*     -> steel blue (HCAL iron absorbers)
            //   anything else  -> light gray (fallback)
            Color_t col      = kGray + 2;
            Char_t  alpha    = 92;   // very transparent so event content shows through
            if      (vname.find("Lead")   != std::string::npos) { col = kGray + 3;   alpha = 88; }
            else if (vname.find("WidePVT") != std::string::npos){ col = kOrange + 1; alpha = 85; }
            else if (vname.find("ThinPS") != std::string::npos) { col = kYellow - 7; alpha = 85; }
            else if (vname.find("Sharp")  != std::string::npos) { col = kAzure - 4;  alpha = 85; }
            else if (vname.find("HPL")    != std::string::npos) { col = kCyan + 2;   alpha = 85; }
            else if (vname.find("HCAL")   != std::string::npos) { col = kGreen + 2;  alpha = 88; }
            else if (vname.find("Iron")   != std::string::npos) { col = kAzure - 2;  alpha = 90; }
            eshape->SetMainColor(col);
            eshape->SetMainTransparency(alpha);
            parent->AddElement(eshape);
        }
        return;  // do not descend
    }
    for (int i = 0; i < node->GetNdaughters(); ++i)
        BuildGeoShapes(node->GetDaughter(i), global, patterns,
                       depth + 1, max_depth, parent);
}

void EventDisplay::OpenEventFile(const std::string& root_path,
                                 const std::string& calo_tree,
                                 const std::string& truth_tree)
{
    fFileCalo = TFile::Open(root_path.c_str(), "READ");
    if (!fFileCalo || fFileCalo->IsZombie()) {
        std::cerr << "[ED] cannot open " << root_path << "\n";
        return;
    }
    fTreeCalo = static_cast<TTree*>(fFileCalo->Get(calo_tree.c_str()));
    fTreeTrue = static_cast<TTree*>(fFileCalo->Get(truth_tree.c_str()));
    if (!fTreeCalo) std::cerr << "[ED] tree '" << calo_tree << "' not found\n";
    if (!fTreeTrue) std::cerr << "[ED] tree '" << truth_tree << "' not found\n";
    if (!fTreeCalo || !fTreeTrue) return;

    fTreeCalo->SetBranchAddress("edep",     &fEdep);
    fTreeCalo->SetBranchAddress("x_global", &fXg);
    fTreeCalo->SetBranchAddress("y_global", &fYg);
    fTreeCalo->SetBranchAddress("z_global", &fZg);

    fTreeTrue->SetBranchAddress("vx_mm",            &fVx_mm);
    fTreeTrue->SetBranchAddress("vy_mm",            &fVy_mm);
    fTreeTrue->SetBranchAddress("vz_mm",            &fVz_mm);
    fTreeTrue->SetBranchAddress("px_MeV",           &fPx_MeV);
    fTreeTrue->SetBranchAddress("py_MeV",           &fPy_MeV);
    fTreeTrue->SetBranchAddress("pz_MeV",           &fPz_MeV);
    fTreeTrue->SetBranchAddress("firstGenAncestor", &fFirstGen);

    std::cout << "[ED] " << fTreeCalo->GetEntries() << " calo events, "
              << fTreeTrue->GetEntries() << " truth events\n";
}

void EventDisplay::OpenClusterFile(const std::string& cluster_path,
                                   const std::string& tree_name)
{
    fFileClus = TFile::Open(cluster_path.c_str(), "READ");
    if (!fFileClus || fFileClus->IsZombie()) {
        std::cerr << "[ED] cannot open " << cluster_path << "\n";
        return;
    }
    fTreeClus = static_cast<TTree*>(fFileClus->Get(tree_name.c_str()));
    if (!fTreeClus) {
        std::cerr << "[ED] tree '" << tree_name << "' not found\n";
        return;
    }
    fTreeClus->SetBranchAddress("hit_ori",     &fHitOri);
    fTreeClus->SetBranchAddress("hit_pos",     &fHitPos);
    fTreeClus->SetBranchAddress("hit_mip",     &fHitMip);
    fTreeClus->SetBranchAddress("hit_z",       &fHitZ);
    fTreeClus->SetBranchAddress("hit_weight",  &fHitWeight);
    fTreeClus->SetBranchAddress("hit_cluster", &fHitClust);
    fTreeClus->SetBranchAddress("cl_id",       &fClId);
    fTreeClus->SetBranchAddress("cl_ori",      &fClOri);
    fTreeClus->SetBranchAddress("cl_mip",      &fClMip);
    fTreeClus->SetBranchAddress("match_h_cl",  &fMatchHcl);
    fTreeClus->SetBranchAddress("match_v_cl",  &fMatchVcl);
    fTreeClus->SetBranchAddress("match_score", &fMatchScore);
    std::cout << "[ED] cluster tree: " << fTreeClus->GetEntries() << " entries\n";
}

void EventDisplay::AutoCenterOnFirstEvent()
{
    // Prefer calo_events hits (the rich edep cloud); fall back to
    // cluster_tree hit_z if calo_events isn't available.
    double xmin = 1e30, xmax = -1e30, ymin = 1e30, ymax = -1e30, zmin = 1e30, zmax = -1e30;
    bool got = false;
    if (fTreeCalo && fTreeCalo->GetEntries() > 0) {
        fTreeCalo->GetEntry(0);
        if (fEdep && !fEdep->empty()) {
            for (size_t k = 0; k < fEdep->size(); ++k) {
                double x = (*fXg)[k] * fHitScale;
                double y = (*fYg)[k] * fHitScale;
                double z = (*fZg)[k] * fHitScale;
                xmin = std::min(xmin, x); xmax = std::max(xmax, x);
                ymin = std::min(ymin, y); ymax = std::max(ymax, y);
                zmin = std::min(zmin, z); zmax = std::max(zmax, z);
            }
            got = true;
        }
    }
    if (!got && fTreeClus && fTreeClus->GetEntries() > 0) {
        fTreeClus->GetEntry(0);
        if (fHitZ && !fHitZ->empty()) {
            for (size_t k = 0; k < fHitZ->size(); ++k) {
                double z = (*fHitZ)[k] * fHitScale;
                zmin = std::min(zmin, z); zmax = std::max(zmax, z);
            }
            xmin = -200.0; xmax = 200.0;   // ad-hoc xy bbox (cm)
            ymin = -200.0; ymax = 200.0;
            got = true;
        }
    }
    if (!got) {
        std::cerr << "[ED] AutoCenterOnFirstEvent: no hits to center on\n";
        return;
    }
    // The scene offset (fOff) is the rigid translation that recenters the
    // view on the detector.  It is computed from the geometry bounding box
    // in AutoCenterGeometry(), which must run after this function.  Here we
    // only report the first-event hit bbox for diagnostics; fOff is left
    // untouched (it will be filled in by AutoCenterGeometry).
    std::cout << "[ED] first-event hit bbox (mm): x[" << xmin << ", " << xmax
              << "] y[" << ymin << ", " << ymax
              << "] z[" << zmin << ", " << zmax << "]\n";
}

// Place the GDML detector geometry into the scene, recentered so the
// detector bounding-box centre sits at the scene origin (0,0,0).
//
// Coordinate handling:
//
//  1) UNIT: TGeo geometry is in centimetres -- both node positions and
//     shape dimensions (TGeoTube/TGeoBBox extents).  The scene is in
//     millimetres, so geometry is scaled by 10 (cm -> mm): node positions
//     are multiplied here, and shape dimensions via SetScale(10,10,10).
//
//  2) FRAME: the geometry is built in the GDML world frame; hits and truth
//     are in the analysis frame.  The analysis-frame origin sits at GDML
//     z = -120827/2 = -60413.5 mm, so  z_analysis = z_gdml + 60413.5 mm.
//
//  3) RECENTRE: ROOT 6.38's web REve auto-fits the camera to the scene
//     bounding box and offers no manual camera/clip control.  Content far
//     from the world origin (the calo is at z ~ 96000 mm) clips and is
//     hard to frame.  So a single rigid translation fOff -- the geometry
//     bbox centre, in analysis-frame mm -- is subtracted from EVERY world
//     coordinate (geometry here, and hits/vertices via toScene()).  This
//     is a pure change of origin: all distances, sizes and angles are
//     preserved exactly.  After it, the detector centre is at the scene
//     origin and the camera frames it by default.
//
// fOff is COMPUTED here (geometry bbox centre) and then used by toScene()
// for everything else, so this must run before any event is drawn.
void EventDisplay::AutoCenterGeometry()
{
    if (fGeoHolder->NumChildren() == 0) return;

    // GDML -> analysis frame: pure z translation, in millimetres.
    const double kFrameShiftMm = 60413.5;
    // cm -> mm scale for geometry (positions AND shape dimensions).
    const double kCmToMm = 10.0;
    // Uniform scene downscale (fHitScale).  Applied to EVERY world quantity
    // -- geometry positions, geometry shape dimensions, and (via toScene)
    // hits / cylinders / vertices.  This is a similarity transform: every
    // distance, size and angle is preserved exactly; only the absolute
    // numbers shrink.  It exists because ROOT 6.38's web REve gives no
    // camera clip-plane control, and a perspective camera cannot resolve a
    // ~58 m scene span (depth-buffer precision collapses -> blank / banded
    // 3D view).  Shrinking the whole scene by fHitScale restores a sane
    // depth range so the 3D perspective view renders.
    const double s = fHitScale;

    // ---- Pass 1: geometry bounding box in SCALED, analysis-frame units --
    double xmin=1e30,xmax=-1e30,ymin=1e30,ymax=-1e30,zmin=1e30,zmax=-1e30;
    for (auto* child : fGeoHolder->RefChildren()) {
        auto* gs = dynamic_cast<REX::REveGeoShape*>(child);
        if (!gs) continue;
        double tx, ty, tz;                       // GDML-frame position, cm
        gs->RefMainTrans().GetPos(tx, ty, tz);
        TGeoBBox* bb = dynamic_cast<TGeoBBox*>(gs->GetShape());
        const double rx = (bb ? bb->GetDX() : 0.0) * kCmToMm * s;
        const double ry = (bb ? bb->GetDY() : 0.0) * kCmToMm * s;
        const double rz = (bb ? bb->GetDZ() : 0.0) * kCmToMm * s;
        // Node position: cm -> mm, frame shift, then uniform scale.
        const double px = (tx * kCmToMm)                 * s;
        const double py = (ty * kCmToMm)                 * s;
        const double pz = (tz * kCmToMm + kFrameShiftMm) * s;
        xmin=std::min(xmin,px-rx); xmax=std::max(xmax,px+rx);
        ymin=std::min(ymin,py-ry); ymax=std::max(ymax,py+ry);
        zmin=std::min(zmin,pz-rz); zmax=std::max(zmax,pz+rz);
    }
    // fOff = detector bbox centre, in SCALED units.  Subtracting it from
    // every (already-scaled) world coordinate puts the detector centre at
    // the origin.  toScene() does mm*fHitScale - fOff, so storing fOff in
    // scaled units keeps it consistent for hits/cylinders/vertices too.
    fOff[0] = 0.5 * (xmin + xmax);
    fOff[1] = 0.5 * (ymin + ymax);
    fOff[2] = 0.5 * (zmin + zmax);
    fOffValid = true;

    // ---- Pass 2: place each geometry node, scaled and recentered -------
    for (auto* child : fGeoHolder->RefChildren()) {
        auto* gs = dynamic_cast<REX::REveGeoShape*>(child);
        if (!gs) continue;
        REX::REveTrans& tr = gs->RefMainTrans();
        double tx, ty, tz;                       // GDML-frame position, cm
        tr.GetPos(tx, ty, tz);
        // Shape dimensions: cm -> mm -> scaled, all uniform.
        tr.SetScale(kCmToMm * s, kCmToMm * s, kCmToMm * s);
        tr.SetPos((tx * kCmToMm)                 * s - fOff[0],
                  (ty * kCmToMm)                 * s - fOff[1],
                  (tz * kCmToMm + kFrameShiftMm) * s - fOff[2]);
    }
    fGeoOffValid = true;
    std::cout << "[ED] geometry recentered on detector; scene offset fOff "
                 "(analysis-frame mm) = ("
              << fOff[0] << ", " << fOff[1] << ", " << fOff[2] << ")\n"
              << "     detector bbox (analysis mm): x[" << xmin << ", " << xmax
              << "] y[" << ymin << ", " << ymax
              << "] z[" << zmin << ", " << zmax << "]\n";
}

// =============================================================================
//  Scene-placement diagnostics
//
//  Prints, for the requested event, the bounding boxes / positions of:
//    * the GDML geometry (as actually placed in the geometry holder),
//    * the calo edep hit cloud (calo_events),
//    * the cluster hits (cluster_tree),
//    * the truth vertex and the reco POCA vertex.
//  The scene is scaled by fHitScale and recentered on the detector (rigid
//  translation by fOff), so "scene" values are small and centred near 0
//  for the detector.  Each line also reports the TRUE analysis-frame mm
//  value -- recovered as (scene + fOff) / fHitScale -- so the physical
//  coordinates remain directly inspectable.
// =============================================================================
void EventDisplay::DumpScenePlacement(Long64_t i)
{
    auto toScene = [&](double xm, double ym, double zm) {
        return Vec3{ xm * fHitScale - fOff[0],
                     ym * fHitScale - fOff[1],
                     zm * fHitScale - fOff[2] };
    };
    // Recover a true analysis-frame mm coordinate from a scaled scene one.
    auto toMm = [&](double scene, int axis) {
        return (scene + fOff[axis]) / fHitScale;
    };
    auto printBox = [&](const char* tag,
                        double x0, double x1, double y0, double y1,
                        double z0, double z1) {
        // z0,z1 are scaled scene coords; report scene centre and true mm.
        const double zc_scene = 0.5*(z0+z1);
        std::cout << "  [place] " << tag
                  << "  scene x[" << x0 << ", " << x1 << "]"
                  << " y[" << y0 << ", " << y1 << "]"
                  << " z[" << z0 << ", " << z1 << "]"
                  << "  (scene z ctr " << zc_scene
                  << ", true z ctr " << toMm(zc_scene, 2) << " mm)\n";
    };

    std::cout << "[ED] ---- scene placement ----  fOff (analysis mm) = ("
              << fOff[0] << ", " << fOff[1] << ", " << fOff[2] << ")\n";

    // --- Geometry: read the placed transforms + shape half-extents -----
    {
        double x0=1e30,x1=-1e30,y0=1e30,y1=-1e30,z0=1e30,z1=-1e30;
        bool any = false;
        for (auto* child : fGeoHolder->RefChildren()) {
            auto* gs = dynamic_cast<REX::REveGeoShape*>(child);
            if (!gs) continue;
            double tx, ty, tz;
            gs->RefMainTrans().GetPos(tx, ty, tz);   // placed scaled scene
            TGeoBBox* bb = dynamic_cast<TGeoBBox*>(gs->GetShape());
            // Shape half-extents are in TGeo cm; the geometry transform
            // carries a uniform scale of (cm->mm) * fHitScale = 10*fHitScale,
            // so the half-extents scale by the same factor.
            const double s = 10.0 * fHitScale;
            const double rx = bb ? bb->GetDX()*s : 0.0;
            const double ry = bb ? bb->GetDY()*s : 0.0;
            const double rz = bb ? bb->GetDZ()*s : 0.0;
            x0=std::min(x0,tx-rx); x1=std::max(x1,tx+rx);
            y0=std::min(y0,ty-ry); y1=std::max(y1,ty+ry);
            z0=std::min(z0,tz-rz); z1=std::max(z1,tz+rz);
            any = true;
        }
        if (any) printBox("GEOMETRY    ", x0,x1,y0,y1,z0,z1);
        else     std::cout << "  [place] GEOMETRY     : no shapes\n";
    }

    // --- Calo edep hits (calo_events) ---------------------------------
    if (fTreeCalo && i < fTreeCalo->GetEntries()) {
        fTreeCalo->GetEntry(i);
        if (fEdep && !fEdep->empty()) {
            double x0=1e30,x1=-1e30,y0=1e30,y1=-1e30,z0=1e30,z1=-1e30;
            for (size_t k = 0; k < fEdep->size(); ++k) {
                Vec3 s = toScene((*fXg)[k], (*fYg)[k], (*fZg)[k]);
                x0=std::min(x0,s.x); x1=std::max(x1,s.x);
                y0=std::min(y0,s.y); y1=std::max(y1,s.y);
                z0=std::min(z0,s.z); z1=std::max(z1,s.z);
            }
            printBox("CALO HITS   ", x0,x1,y0,y1,z0,z1);
        } else {
            std::cout << "  [place] CALO HITS    : no edep hits this event\n";
        }
    }

    // --- Cluster hits (cluster_tree) ----------------------------------
    if (fTreeClus && i < fTreeClus->GetEntries()) {
        fTreeClus->GetEntry(i);
        if (fHitZ && !fHitZ->empty()) {
            // hit_pos carries only one transverse coordinate per hit
            // (x for V hits, y for H hits); we still bracket z fully.
            double z0=1e30,z1=-1e30, p0=1e30,p1=-1e30;
            for (size_t k = 0; k < fHitZ->size(); ++k) {
                double zc = (*fHitZ)[k]   * fHitScale - fOff[2];
                double pc = (*fHitPos)[k] * fHitScale;   // transverse, no fixed offset
                z0=std::min(z0,zc); z1=std::max(z1,zc);
                p0=std::min(p0,pc); p1=std::max(p1,pc);
            }
            std::cout << "  [place] CLUSTER HITS "
                      << "  z[" << z0 << ", " << z1 << "]"
                      << "  (center z = " << 0.5*(z0+z1) << ")"
                      << "  transverse pos[" << p0 << ", " << p1 << "]\n";
        } else {
            std::cout << "  [place] CLUSTER HITS : no cluster hits this event\n";
        }
    }

    // --- Truth + reco vertices ----------------------------------------
    if (fTreeTrue && i < fTreeTrue->GetEntries()) fTreeTrue->GetEntry(i);
    if (fTreeClus && i < fTreeClus->GetEntries()) fTreeClus->GetEntry(i);
    Vec3 truthVtx; Vec3 truthDir[2];
    Vec3 recoVtx;  PhotonReco phs[2];
    const bool tOk = TruthFromCurrentEvent(truthVtx, truthDir);
    const bool rOk = ReconstructCurrentEvent(phs, recoVtx);
    if (tOk) {
        Vec3 s = toScene(truthVtx.x, truthVtx.y, truthVtx.z);
        std::cout << "  [place] TRUTH VERTEX "
                  << "  mm(" << truthVtx.x << ", " << truthVtx.y << ", " << truthVtx.z << ")"
                  << "  ->  scene(" << s.x << ", " << s.y << ", " << s.z << ")\n";
    } else {
        std::cout << "  [place] TRUTH VERTEX : unavailable this event\n";
    }
    if (rOk) {
        Vec3 s = toScene(recoVtx.x, recoVtx.y, recoVtx.z);
        std::cout << "  [place] RECO  POCA   "
                  << "  mm(" << recoVtx.x << ", " << recoVtx.y << ", " << recoVtx.z << ")"
                  << "  ->  scene(" << s.x << ", " << s.y << ", " << s.z << ")\n";
    } else {
        std::cout << "  [place] RECO  POCA   : unavailable this event\n";
    }
    std::cout << "[ED] --------------------------------------------------\n";
}

// =============================================================================
//  Truth + reco extraction
// =============================================================================
bool EventDisplay::TruthFromCurrentEvent(Vec3& vtx, Vec3 dir[2]) const
{
    if (!fFirstGen) return false;
    int pi0Idx = -1, pIdx[2] = { -1, -1 };
    for (size_t i = 0; i < fFirstGen->size(); ++i) {
        int fga = (*fFirstGen)[i];
        if (fga == 1 && pi0Idx < 0) pi0Idx = (int)i;
        if (fga == 2 && pIdx[0] < 0) pIdx[0] = (int)i;
        if (fga == 3 && pIdx[1] < 0) pIdx[1] = (int)i;
    }
    if (pIdx[0] < 0 || pIdx[1] < 0) return false;
    const int vtxIdx = (pi0Idx >= 0) ? pi0Idx : pIdx[0];
    vtx = { (*fVx_mm)[vtxIdx], (*fVy_mm)[vtxIdx], (*fVz_mm)[vtxIdx] };
    for (int p = 0; p < 2; ++p) {
        double px = (*fPx_MeV)[pIdx[p]];
        double py = (*fPy_MeV)[pIdx[p]];
        double pz = (*fPz_MeV)[pIdx[p]];
        double m  = std::sqrt(px*px + py*py + pz*pz);
        if (m <= 0) return false;
        dir[p] = { px/m, py/m, pz/m };
    }
    return true;
}

bool EventDisplay::ReconstructCurrentEvent(PhotonReco phs[2], Vec3& pocaVtx) const
{
    if (!fMatchScore || fMatchScore->size() < 2) return false;

    // Pick top-2 matches by lowest score
    std::vector<std::pair<double, int>> rank;
    for (size_t m = 0; m < fMatchScore->size(); ++m)
        rank.push_back({ (*fMatchScore)[m], (int)m });
    std::sort(rank.begin(), rank.end(),
              [](auto& a, auto& b){ return a.first < b.first; });
    int mi[2] = { rank[0].second, rank[1].second };

    const double rMaxMm    = 50.0;
    const int    maxIter   = 10;
    const double convTolMm = 1.0;

    for (int p = 0; p < 2; ++p) {
        int hId = (*fMatchHcl)[mi[p]];
        int vId = (*fMatchVcl)[mi[p]];
        phs[p].hClId = hId;
        phs[p].vClId = vId;
        std::vector<double> Zh, Yh, Wh, Zv, Xv, Wv;
        for (size_t h = 0; h < fHitClust->size(); ++h) {
            int hc = (*fHitClust)[h];
            double w = (*fHitMip)[h] * (*fHitWeight)[h];
            if (hc == hId) {
                Zh.push_back((*fHitZ)[h]);
                Yh.push_back((*fHitPos)[h]);
                Wh.push_back(w);
            } else if (hc == vId) {
                Zv.push_back((*fHitZ)[h]);
                Xv.push_back((*fHitPos)[h]);
                Wv.push_back(w);
            }
        }
        if (Zh.size() < 3 || Zv.size() < 3) return false;

        LineFit2D Fx = iterateHaloCut2D(Zv, Xv, Wv, rMaxMm, maxIter, convTolMm);
        LineFit2D Fy = iterateHaloCut2D(Zh, Yh, Wh, rMaxMm, maxIter, convTolMm);
        phs[p].line = combineFits(Fx, Fy);
        if (!phs[p].line.ok) return false;

        // MIP-weighted 3D centroid: use H cluster hits for (y, z) and
        // V cluster hits for (x, z); average z across the two cluster axes.
        double sumW_h = 0, sumWy = 0, sumWz_h = 0;
        for (size_t i = 0; i < Zh.size(); ++i) {
            sumW_h += Wh[i]; sumWy += Wh[i] * Yh[i]; sumWz_h += Wh[i] * Zh[i];
        }
        double sumW_v = 0, sumWx = 0, sumWz_v = 0;
        for (size_t i = 0; i < Zv.size(); ++i) {
            sumW_v += Wv[i]; sumWx += Wv[i] * Xv[i]; sumWz_v += Wv[i] * Zv[i];
        }
        if (sumW_h <= 0 || sumW_v <= 0) return false;
        const double yc = sumWy / sumW_h;
        const double xc = sumWx / sumW_v;
        const double zc = 0.5 * (sumWz_h / sumW_h + sumWz_v / sumW_v);
        phs[p].centroid = { xc, yc, zc };

        // Cluster z extent (half-length): use the min/max z of all hits.
        double zmin = 1e30, zmax = -1e30;
        for (double z : Zh) { zmin = std::min(zmin, z); zmax = std::max(zmax, z); }
        for (double z : Zv) { zmin = std::min(zmin, z); zmax = std::max(zmax, z); }
        phs[p].z_extent_half = 0.5 * (zmax - zmin);

        // Transverse RMS about the fitted line.
        // Project each hit's residual onto the plane perpendicular to the line.
        // For an H hit only y is measured: transverse y-residual = y - y_pred(z).
        // For a V hit only x is measured: transverse x-residual = x - x_pred(z).
        // Treat them separately; combine to a 2D "spread" by sqrt(rms_x^2 + rms_y^2)/sqrt(2).
        double sumWyy = 0, sumWxx = 0;
        double sw_y = 0, sw_x = 0;
        for (size_t i = 0; i < Zh.size(); ++i) {
            double y_pred = Fy.intercept + Fy.slope * (Zh[i] - Fy.zRef);
            double dy = Yh[i] - y_pred;
            sumWyy += Wh[i] * dy * dy;
            sw_y   += Wh[i];
        }
        for (size_t i = 0; i < Zv.size(); ++i) {
            double x_pred = Fx.intercept + Fx.slope * (Zv[i] - Fx.zRef);
            double dx = Xv[i] - x_pred;
            sumWxx += Wv[i] * dx * dx;
            sw_x   += Wv[i];
        }
        double rms_y = sw_y > 0 ? std::sqrt(sumWyy / sw_y) : 0.0;
        double rms_x = sw_x > 0 ? std::sqrt(sumWxx / sw_x) : 0.0;
        // 2D combined transverse rms (averaging quadratically across axes)
        phs[p].r_rms = std::sqrt(0.5 * (rms_x*rms_x + rms_y*rms_y));
        phs[p].ok = true;
    }

    auto pres = poca(phs[0].line, phs[1].line);
    if (!pres.ok) return false;
    pocaVtx = pres.vtx;
    return true;
}

// =============================================================================
//  Best-event picker
// =============================================================================
Long64_t EventDisplay::PickBestEvent()
{
    if (!fTreeClus || !fTreeTrue) {
        std::cerr << "[ED] cluster/truth trees missing -- cannot pick best event\n";
        return 0;
    }
    const Long64_t N = std::min(fTreeClus->GetEntries(), fTreeTrue->GetEntries());
    Long64_t bestEv = 0;
    double bestRes = std::numeric_limits<double>::infinity();
    for (Long64_t i = 0; i < N; ++i) {
        fTreeClus->GetEntry(i);
        fTreeTrue->GetEntry(i);
        Vec3 truthVtx; Vec3 truthDir[2];
        if (!TruthFromCurrentEvent(truthVtx, truthDir)) continue;
        PhotonReco phs[2]; Vec3 vtx;
        if (!ReconstructCurrentEvent(phs, vtx)) continue;
        double res = std::abs(vtx.z - truthVtx.z);
        if (res < bestRes) { bestRes = res; bestEv = i; }
        if ((i + 1) % 200 == 0)
            std::cout << "[ED] scanned " << (i+1) << "/" << N
                      << "  best so far: ev " << bestEv
                      << "  |Dz|=" << bestRes << " mm\n";
    }
    std::cout << "[ED] best event: " << bestEv
              << "  |z_reco - z_truth| = " << bestRes << " mm\n";
    return bestEv;
}

// =============================================================================
//  Per-event display
// =============================================================================
void EventDisplay::GotoEvent(Long64_t i)
{
    // Choose a representative tree for bounds-checking; prefer cluster
    // tree (the one carrying reco), fall back to calo, then truth.
    TTree* refTree = fTreeClus ? fTreeClus : (fTreeCalo ? fTreeCalo : fTreeTrue);
    if (!refTree) {
        std::cerr << "[ED] no trees loaded\n";
        return;
    }
    if (i < 0 || i >= refTree->GetEntries()) {
        std::cerr << "[ED] event " << i << " out of range\n";
        return;
    }
    if (fTreeCalo && i < fTreeCalo->GetEntries()) fTreeCalo->GetEntry(i);
    if (fTreeTrue && i < fTreeTrue->GetEntries()) fTreeTrue->GetEntry(i);
    if (fTreeClus && i < fTreeClus->GetEntries()) fTreeClus->GetEntry(i);

    REX::REveManager::ChangeGuard guard;
    fEventHolder->DestroyElements();

    // ------------------------------------------------------------------
    // 1) Backdrop: edep point cloud (5 log-binned colors).  Skipped if
    //    the calo_events tree didn't load or this event has no hits.
    // ------------------------------------------------------------------
    if (fEdep && !fEdep->empty()) {
        const double emin = *std::min_element(fEdep->begin(), fEdep->end());
        const double emax = std::max(emin + 1e-6,
                                     *std::max_element(fEdep->begin(), fEdep->end()));
        constexpr int kNBins = 5;
        const Color_t bin_colors[kNBins] = {
            kAzure + 1, kCyan + 1, kGreen + 1, kOrange + 7, kRed + 1
        };
        const Float_t bin_sizes[kNBins] = { 4.5f, 5.5f, 6.5f, 7.5f, 8.5f };
        REX::REvePointSet* pbin[kNBins];
        for (int b = 0; b < kNBins; ++b) {
            pbin[b] = new REX::REvePointSet(
                Form("hits_bin_%d", b), Form("Energy bin %d", b));
            pbin[b]->SetMarkerStyle(20);
            pbin[b]->SetMarkerSize(bin_sizes[b]);
            pbin[b]->SetMarkerColor(bin_colors[b]);
        }
        const double log_lo = std::log10(std::max(emin, 1e-9));
        const double log_hi = std::log10(std::max(emax, log_lo * 10));
        const double log_w  = std::max(log_hi - log_lo, 1e-6);
        for (size_t k = 0; k < fEdep->size(); ++k) {
            const float x = (float)((*fXg)[k] * fHitScale - fOff[0]);
            const float y = (float)((*fYg)[k] * fHitScale - fOff[1]);
            const float z = (float)((*fZg)[k] * fHitScale - fOff[2]);
            const double le = std::log10(std::max((*fEdep)[k], 1e-9));
            int b = (int)(((le - log_lo) / log_w) * kNBins);
            if (b < 0) b = 0;
            if (b >= kNBins) b = kNBins - 1;
            pbin[b]->SetNextPoint(x, y, z);
        }
        for (int b = 0; b < kNBins; ++b) fEventHolder->AddElement(pbin[b]);
    } else {
        std::cout << "[ED] event " << i << ": no calo edep hits to draw\n";
    }


    // ------------------------------------------------------------------
    // 2) Truth + reco overlays
    // ------------------------------------------------------------------
    Vec3 truthVtx; Vec3 truthDir[2];
    Vec3 recoVtx;  PhotonReco phs[2];
    const bool tOk = TruthFromCurrentEvent(truthVtx, truthDir);
    const bool rOk = ReconstructCurrentEvent(phs, recoVtx);

    // Helper: convert a mm-space point to scene-space cm point.
    auto toScene = [&](double xm, double ym, double zm) {
        return Vec3{
            xm * fHitScale - fOff[0],
            ym * fHitScale - fOff[1],
            zm * fHitScale - fOff[2]
        };
    };

    // ------------------------------------------------------------------
    // Cylinders along the reco-line directions, centered on each
    // photon's cluster center (the through-point of the fitted line,
    // i.e. the same point the POCA calculation uses as line origin).
    //
    // Implementation:
    //   - One TGeoTube per photon, radius = transverse MIP-weighted RMS,
    //     half-length = half of the cluster z extent (along the line).
    //   - Wrapped in a REveGeoShape with a TGeoCombiTrans transform that
    //     rotates the tube's z-axis to align with the reco line direction
    //     and translates to the cluster center.
    //   - Also add a small bright marker at the cluster center so any
    //     positioning issue is immediately visible (the dot tells us
    //     where the cylinder SHOULD be; if the dot is visible but the
    //     cylinder isn't, that's a rendering issue, not a placement one).
    // ------------------------------------------------------------------
    if (rOk) {
        const Color_t cylCol[2] = { kCyan + 1, kMagenta + 1 };

        // Full calo z extent from ALL cluster hits in this event.  Used as
        // the cylinder length so each photon's tube spans the entire
        // detector volume rather than just its own cluster's reach.
        double caloZmin_mm = 1e30, caloZmax_mm = -1e30;
        if (fHitZ) {
            for (size_t h = 0; h < fHitZ->size(); ++h) {
                const double z = (*fHitZ)[h];
                if (z < caloZmin_mm) caloZmin_mm = z;
                if (z > caloZmax_mm) caloZmax_mm = z;
            }
        }
        const double caloHalfLen_mm = (caloZmax_mm > caloZmin_mm)
                                       ? 0.5 * (caloZmax_mm - caloZmin_mm)
                                       : 750.0;   // fallback ~1.5 m

        for (int p = 0; p < 2; ++p) {
            if (!phs[p].ok) continue;

            // Cluster center: project the line through-point onto the
            // center of the calo z range (so the cylinder spans the calo
            // symmetrically when oriented along the reco line direction).
            const double caloZmid_mm = 0.5 * (caloZmin_mm + caloZmax_mm);
            const Vec3   d  = phs[p].line.d;
            const double inv_dz = (std::abs(d.z) > 1e-9) ? 1.0 / d.z : 0.0;
            // Project the line origin along the line direction to caloZmid:
            //   c = line.p + t * line.d  where t scales so c.z = caloZmid
            const double t_mm = (caloZmid_mm - phs[p].line.p.z) * inv_dz;
            const Vec3 c_mm{
                phs[p].line.p.x + d.x * t_mm,
                phs[p].line.p.y + d.y * t_mm,
                caloZmid_mm
            };

            // Inflate radius for visibility (physical r_rms is often a few
            // mm; 4x makes it visible alongside the calo at ~50 mm pitch).
            const double hl_mm = caloHalfLen_mm;
            const double r_mm  = std::max(80.0, 4.0 * phs[p].r_rms);

            // Convert to scene-cm.
            const double hl_cm = hl_mm * fHitScale;
            const double r_cm  = r_mm  * fHitScale;
            const Vec3   c_cm  = toScene(c_mm.x, c_mm.y, c_mm.z);

            // ---- Build rotation that maps z_hat -> d_hat -----------------
            // Using axis-angle form: rotation axis u = (z_hat x d_hat),
            // angle = acos(d.z).  If d is already (very near) the z-axis,
            // skip the rotation -- the small-angle case dominates anyway
            // for forward-going photons.
            double rot[9] = { 1, 0, 0,    // row-major 3x3, identity by default
                              0, 1, 0,
                              0, 0, 1 };
            const double dz_clamped = std::max(-1.0, std::min(1.0, d.z));
            const double ang_rad    = std::acos(dz_clamped);
            // Axis u = z x d = (-d.y, d.x, 0).
            const double ux = -d.y;
            const double uy =  d.x;
            const double uz =  0.0;
            const double u_norm = std::sqrt(ux*ux + uy*uy + uz*uz);
            if (u_norm > 1e-9 && ang_rad > 1e-6) {
                const double inv = 1.0 / u_norm;
                const double Ux = ux * inv;
                const double Uy = uy * inv;
                const double Uz = uz * inv;
                const double s  = std::sin(ang_rad);
                const double cc = std::cos(ang_rad);
                const double C  = 1.0 - cc;
                // Rodrigues' rotation matrix, row-major.
                rot[0] = cc + Ux*Ux*C;
                rot[1] = Ux*Uy*C - Uz*s;
                rot[2] = Ux*Uz*C + Uy*s;
                rot[3] = Uy*Ux*C + Uz*s;
                rot[4] = cc + Uy*Uy*C;
                rot[5] = Uy*Uz*C - Ux*s;
                rot[6] = Uz*Ux*C - Uy*s;
                rot[7] = Uz*Uy*C + Ux*s;
                rot[8] = cc + Uz*Uz*C;
            }

            // Sanity check the rotation by applying it to z_hat and
            // confirming we get d back.  Print result for diagnostics.
            const double rz_x = rot[2];
            const double rz_y = rot[5];
            const double rz_z = rot[8];
            std::cout << "  photon " << p
                      << "  cylinder center (scene cm) = ("
                      << c_cm.x << ", " << c_cm.y << ", " << c_cm.z << ")"
                      << "  half-len=" << hl_cm << " cm  r=" << r_cm << " cm\n"
                      << "    line d = (" << d.x << ", " << d.y << ", " << d.z << ")"
                      << "    R*z_hat = (" << rz_x << ", " << rz_y << ", " << rz_z << ")\n";

            // ---- Build TGeoCombiTrans (translation + rotation) -----------
            TGeoRotation rotation;
            rotation.SetMatrix(rot);
            TGeoTranslation translation(c_cm.x, c_cm.y, c_cm.z);
            // Combine: TGeoCombiTrans constructor takes translation then
            // rotation, applying rotation first, then translation.
            TGeoCombiTrans combi(translation, rotation);

            // ---- Build the shape (TGeoTube) and wrap it ------------------
            // Use a unique name per shape; REveGeoShape needs a TGeoShape*.
            // The TGeoTube is owned by REveGeoShape after SetShape().
            auto* tube = new TGeoTube(0.0, r_cm, hl_cm);

            auto* shape = new REX::REveGeoShape(Form("cluster_cyl_%d", p));
            shape->SetShape(tube);
            shape->SetMainColor(cylCol[p]);
            shape->SetMainTransparency(20);
            shape->RefMainTrans().SetFrom(combi);
            fEventHolder->AddElement(shape);

            // ---- Line from cluster center to the reco (POCA) vertex ------
            if (rOk) {
                auto* ls = new REX::REveStraightLineSet(
                    Form("reco_line_%d", p));
                const Vec3 v = toScene(recoVtx.x, recoVtx.y, recoVtx.z);
                ls->AddLine(c_cm.x, c_cm.y, c_cm.z, v.x, v.y, v.z);
                ls->SetLineColor(cylCol[p]);
                ls->SetLineWidth(4);
                fEventHolder->AddElement(ls);
            }
        }

        // ---- Vertex markers ----------------------------------------------
        // Truth vertex: green sphere.  Reco vertex (POCA): yellow sphere.
        // No lines emanate from the truth vertex per request -- only the
        // marker.  Both vertices use REveGeoShape + TGeoSphere because we
        // know that pipeline works in this REve build.
        auto makeVtxSphere = [&](const Vec3& v_mm, Color_t col, const char* name) {
            const Vec3 v_scene = toScene(v_mm.x, v_mm.y, v_mm.z);
            // Radius: a physical ~400 mm marker, scaled by fHitScale like
            // every other world quantity so it stays proportionate in the
            // downscaled scene (uniform scale -> no distortion).
            const double sphR = 400.0 * fHitScale;
            auto* sph = new TGeoSphere(0.0, sphR);
            auto* gs  = new REX::REveGeoShape(name);
            gs->SetShape(sph);
            gs->SetMainColor(col);
            gs->SetMainTransparency(0);   // fully opaque
            TGeoTranslation tr(v_scene.x, v_scene.y, v_scene.z);
            gs->RefMainTrans().SetFrom(tr);
            fEventHolder->AddElement(gs);
        };
        if (tOk) makeVtxSphere(truthVtx, kGreen + 2,   "truth_vertex");
        if (rOk) makeVtxSphere(recoVtx,  kYellow,      "reco_vertex_POCA");
    }


    std::cout << "[ED] event " << i << ": "
              << fEdep->size() << " calo hits"
              << (rOk ? "; reco ok" : "; reco FAILED")
              << (tOk ? "; truth ok" : "; truth FAILED") << "\n";
    if (rOk && tOk) {
        std::cout << "  truth vtx (mm): (" << truthVtx.x << ", " << truthVtx.y
                  << ", " << truthVtx.z << ")\n";
        std::cout << "  reco POCA (mm): (" << recoVtx.x << ", " << recoVtx.y
                  << ", " << recoVtx.z << ")\n";
        std::cout << "  Dz = " << (recoVtx.z - truthVtx.z) << " mm\n";
        // Scene-space positions (cm) -- useful for diagnosing camera framing
        Vec3 tv = toScene(truthVtx.x, truthVtx.y, truthVtx.z);
        Vec3 rv = toScene(recoVtx.x,  recoVtx.y,  recoVtx.z);
        std::cout << "  truth vtx (scene cm): (" << tv.x << ", " << tv.y
                  << ", " << tv.z << ")\n";
        std::cout << "  reco POCA (scene cm): (" << rv.x << ", " << rv.y
                  << ", " << rv.z << ")\n";
        for (int p = 0; p < 2; ++p) {
            Vec3 c_scene = toScene(phs[p].centroid.x,
                                   phs[p].centroid.y,
                                   phs[p].centroid.z);
            std::cout << "  photon " << p
                      << "  cluster centroid (scene cm): ("
                      << c_scene.x << ", " << c_scene.y << ", "
                      << c_scene.z << ")"
                      << "  r_rms=" << phs[p].r_rms << " mm"
                      << "  half-len=" << phs[p].z_extent_half << " mm\n";
        }
    }
}

// =============================================================================
//  main()
// =============================================================================
namespace {
    EventDisplay* g_display = nullptr;
}

int main(int argc, char** argv)
{
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <detector.gdml> <truth+calo.root> <cluster.root>"
                  << " [--event N] [--scale F]"
                  << " [--tree-truth calo_events]"
                  << " [--tree-particles truth_particles]"
                  << " [--tree-cluster cluster_tree]\n";
        return 1;
    }
    // Capture file path args BEFORE TApplication is constructed --
    // TApplication's constructor takes &argc, &argv and (in some ROOT
    // versions) clobbers entries it doesn't recognise, which is why
    // these would otherwise come out empty.
    const std::string gdml_path     = argv[1];
    const std::string file_path_a   = argv[2];
    const std::string file_path_b   = argv[3];

    std::string calo_tree      = "calo_events";
    std::string particles_tree = "truth_particles";
    std::string cluster_tree   = "cluster_tree";
    Long64_t userEvent = -1;
    double userScale = -1.0;   // <= 0 means use default
    for (int i = 4; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--event" && i + 1 < argc) {
            userEvent = std::atoll(argv[++i]);
        } else if (a == "--scale" && i + 1 < argc) {
            userScale = std::atof(argv[++i]);
        } else if (a == "--tree-truth" && i + 1 < argc) {
            calo_tree = argv[++i];
        } else if (a == "--tree-particles" && i + 1 < argc) {
            particles_tree = argv[++i];
        } else if (a == "--tree-cluster" && i + 1 < argc) {
            cluster_tree = argv[++i];
        }
    }

    TApplication app("eve_app", &argc, argv);

    EventDisplay ed;
    g_display = &ed;
    ed.Init();
    if (userScale > 0) {
        // fHitScale is now applied uniformly to geometry positions, geometry
        // shape dimensions AND hits/cylinders/vertices, so overriding it
        // simply rescales the whole scene -- no desync.  Useful if the
        // default 0.01 needs tuning for the 3D perspective view.
        ed.SetHitScale(userScale);
        std::cout << "[ED] scene scale overridden: fHitScale = "
                  << userScale << "\n";
    }

    // Geometry filters: same as original.
    const std::vector<std::string> include = {
        R"(^ECAL_GL\d+_)",
        R"(^HCAL_GL\d+_)",
        R"(^IronPlate)",
    };
    ed.LoadGeometry(gdml_path, include);
    // AutoCenterGeometry() is deferred until after AutoCenterOnFirstEvent()
    // has run: the geometry is now placed using the shared hit offset
    // (fOff), which is only known once the first event's calo hits are read.

    // Auto-detect which of the two ROOT file arguments holds the
    // (calo_events + truth_particles) trees vs the cluster_tree.  This
    // way the user can pass them in either order on the command line.
    auto fileHasTree = [](const std::string& path, const std::string& tree) {
        TFile* f = TFile::Open(path.c_str(), "READ");
        if (!f || f->IsZombie()) return false;
        TTree* t = static_cast<TTree*>(f->Get(tree.c_str()));
        const bool has = (t != nullptr);
        f->Close();
        delete f;
        return has;
    };

    std::string truth_path   = file_path_a;
    std::string cluster_path = file_path_b;
    const bool a2_has_cluster = fileHasTree(truth_path,   cluster_tree);
    const bool a3_has_cluster = fileHasTree(cluster_path, cluster_tree);
    if (a2_has_cluster && !a3_has_cluster) {
        std::cout << "[ED] file order looks reversed; swapping.\n";
        std::swap(truth_path, cluster_path);
    } else if (!a2_has_cluster && !a3_has_cluster) {
        std::cerr << "[ED] neither file contains tree '" << cluster_tree
                  << "' -- check the input paths\n";
    }

    ed.OpenEventFile(truth_path, calo_tree, particles_tree);
    ed.OpenClusterFile(cluster_path, cluster_tree);
    ed.AutoCenterOnFirstEvent();
    // AutoCenterGeometry() places the detector AND computes fOff -- the
    // rigid translation that recenters the whole scene on the detector.
    // It must run before any event is drawn (toScene() uses fOff).
    ed.AutoCenterGeometry();

    // Pick the best event (smallest |z_reco - z_truth|) unless overridden.
    Long64_t evToShow = 0;
    if (userEvent >= 0) {
        evToShow = userEvent;
        std::cout << "[ED] using user-specified event " << evToShow << "\n";
    } else if (ed.NumEvents() > 0) {
        evToShow = ed.PickBestEvent();
    }
    if (ed.NumEvents() > 0) ed.GotoEvent(evToShow);

    // Print where geometry / hits / vertex actually land in scene space,
    // so any remaining frame misalignment is directly visible.
    if (ed.NumEvents() > 0) ed.DumpScenePlacement(evToShow);

    // Add a viewer whose camera auto-frames the detector geometry.
    ed.SetupViewers();

    ed.Manager()->Show();
    app.Run();
    return 0;
}
