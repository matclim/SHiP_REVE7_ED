// =============================================================================
//  event_display.cxx
//
//  REve-based event display for a sampling calorimeter:
//    * Loads a GDML geometry, keeping ONLY volumes whose names match user-
//      provided regular expressions.  For the example detector this means
//      the 20 ECAL layer envelopes (Lead + WidePVT + ThinPS) and skips the
//      thousands of HPL hodoscope fibers entirely.  This is the trick that
//      keeps the display fast: a few dozen shapes upload to the GPU once
//      and persist for the session.
//    * Reads a TTree of calorimeter hits (per-event std::vector<double>
//      branches: edep, x_global, y_global, z_global) and draws them as a
//      single REveBoxSet -- one GPU-instanced collection per event,
//      coloured by deposited energy via REveRGBAPalette.
//
//  Build with the accompanying CMakeLists.txt against an installed ROOT
//  >= 6.30 configured with -Dwebgui=ON -Droot7=ON -Dgdml=ON.
//
//  Run:   ./event_display detector.gdml events.root
//         (then click "next event" / "prev event" buttons in the browser, or
//          use the helper buttons added below)
// =============================================================================

#include <ROOT/REveManager.hxx>
#include <ROOT/REveScene.hxx>
#include <ROOT/REveViewer.hxx>
#include <ROOT/REveElement.hxx>
#include <ROOT/REveGeoShape.hxx>
#include <ROOT/REveBoxSet.hxx>
#include <ROOT/REvePointSet.hxx>
#include <ROOT/REveRGBAPalette.hxx>
#include <ROOT/REveTrans.hxx>

#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TGeoVolume.h>
#include <TGeoMatrix.h>
#include <TGeoShape.h>
#include <TGeoBBox.h>

#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TColor.h>
#include <TString.h>

#include <iostream>
#include <limits>
#include <vector>
#include <regex>
#include <string>
#include <algorithm>
#include <cmath>

namespace REX = ROOT::Experimental;

// -----------------------------------------------------------------------------
//  EventDisplay
// -----------------------------------------------------------------------------
class EventDisplay {
public:
    EventDisplay()  = default;
    ~EventDisplay() = default;

    void Init();

    // Load GDML and keep only volumes whose names match any of the regexes.
    // When a node matches, its shape is captured with its cumulative global
    // transform and we *do not descend further* -- this is what makes the
    // display fast on big detectors.  Use max_depth >= 0 to cap recursion.
    void LoadGeometry(const std::string& gdml_path,
                      const std::vector<std::string>& include_regex,
                      int max_depth = -1);

    void OpenEventFile(const std::string& root_path,
                       const std::string& tree_name = "calo_events");

    Long64_t NumEvents() const { return fTree ? fTree->GetEntries() : 0; }

    // Switch to event i; safe to call from the REve UI.
    void GotoEvent(Long64_t i);
    void NextEvent() { if (fTree) GotoEvent((fCurrent + 1) % NumEvents()); }
    void PrevEvent() { if (fTree) GotoEvent((fCurrent - 1 + NumEvents()) % NumEvents()); }

    // Hand the manager off to the user; caller drives TApplication::Run().
    REX::REveManager* Manager() { return fEve; }

private:
    void BuildGeoShapes(TGeoNode* node, const TGeoHMatrix& parent_mtx,
                        const std::vector<std::regex>& patterns,
                        int depth, int max_depth,
                        REX::REveElement* parent);

    // ---- REve plumbing -----------------------------------------------------
    REX::REveManager*   fEve            = nullptr;
    REX::REveScene*     fGeoScene       = nullptr;   // alias to global scene
    REX::REveScene*     fEventScene     = nullptr;   // alias to event scene
    REX::REveElement*   fGeoHolder      = nullptr;   // parent of all geometry
    REX::REveElement*   fEventHolder    = nullptr;   // parent of per-event data

    // ---- I/O ---------------------------------------------------------------
    TFile*    fFile    = nullptr;
    TTree*    fTree    = nullptr;
    Long64_t  fCurrent = -1;

    // Multiplicative scale applied to every hit position and to the cell
    // drawing size before they're handed to REve.  TGeo internally works
    // in cm; if your hit coordinates come out of Geant4 in mm, leave this
    // at 0.1.  Set to 1.0 if your data is already in cm.
    double    fHitScale = 0.1;

    // Translations subtracted before passing things to REve.  We keep
    // *separate* offsets for geometry and hits so a misalignment between
    // the GDML's world frame and the simulation's hit-coordinate frame
    // doesn't leave one of them off-camera.  Each is computed
    // independently from its own bbox; the result is that both show up
    // overlapping at the scene origin even when they're physically far
    // apart in the world.
    double    fGeoOff[3] = {0.0, 0.0, 0.0};
    double    fHitOff[3] = {0.0, 0.0, 0.0};
    bool      fAutoCenterGeo  = true;
    bool      fAutoCenterHits = true;

    // After AutoCenterGeometry runs, remember the shifted (scene-space)
    // bbox of the detector so the hit auto-centering can align to it.
    double    fGeoSceneMin[3] = {0.0, 0.0, 0.0};
    double    fGeoSceneMax[3] = {0.0, 0.0, 0.0};
    bool      fGeoBboxValid   = false;

    // When both geometry and hits are present, align the hits' z-front
    // (= z_min) to the geometry's z-front instead of centering on the
    // hits' own bbox.  This is physically what you want for a beam-test
    // calorimeter -- the shower starts at the entry face and develops
    // downstream, so its centroid is in front of the detector centroid.
    // x and y are still centered on each bbox's centre.
    bool      fAlignHitsZFront = true;

public:
    void SetHitScale(double s) { fHitScale = s; }
    void SetAutoCenter(bool b) { fAutoCenterGeo = b; fAutoCenterHits = b; }
    void SetAlignHitsZFront(bool b) { fAlignHitsZFront = b; }
    // Pin the geometry offset; auto-centering for geometry is then off.
    void SetGeoOffset(double x, double y, double z) {
        fGeoOff[0] = x; fGeoOff[1] = y; fGeoOff[2] = z;
        fAutoCenterGeo = false;
    }
    // Pin the hit offset; auto-centering for hits is then off.
    void SetHitOffset(double x, double y, double z) {
        fHitOff[0] = x; fHitOff[1] = y; fHitOff[2] = z;
        fAutoCenterHits = false;
    }
    // Pin both at once (back-compat with previous SetSceneOffset name).
    void SetSceneOffset(double x, double y, double z) {
        SetGeoOffset(x, y, z);
        SetHitOffset(x, y, z);
    }
    // Compute geometry offset from loaded shapes' bbox and shift them.
    void AutoCenterGeometry();
    // Compute hit offset from first event's hits' bbox.
    void AutoCenterOnFirstEvent();
private:

    // ---- Branch buffers (must match calo_events tree) ----------------------
    std::vector<double>* fEdep    = nullptr;
    std::vector<double>* fXg      = nullptr;
    std::vector<double>* fYg      = nullptr;
    std::vector<double>* fZg      = nullptr;
    std::vector<int>*    fType    = nullptr;
    std::vector<int>*    fSection = nullptr;
    std::vector<int>*    fLayer   = nullptr;
    std::vector<int>*    fHcal    = nullptr;
    std::vector<int>*    fHexant  = nullptr;
};

// -----------------------------------------------------------------------------
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

// -----------------------------------------------------------------------------
void EventDisplay::LoadGeometry(const std::string& gdml_path,
                                const std::vector<std::string>& include_regex,
                                int max_depth)
{
    TGeoManager* geo = TGeoManager::Import(gdml_path.c_str());
    if (!geo) {
        std::cerr << "[EventDisplay] failed to import GDML: "
                  << gdml_path << "\n";
        return;
    }

    std::vector<std::regex> patterns;
    patterns.reserve(include_regex.size());
    for (const auto& s : include_regex) {
        try { patterns.emplace_back(s); }
        catch (const std::regex_error& e) {
            std::cerr << "[EventDisplay] bad regex '" << s
                      << "': " << e.what() << "\n";
        }
    }

    {
        REX::REveManager::ChangeGuard guard;
        TGeoHMatrix identity;
        BuildGeoShapes(geo->GetTopNode(), identity, patterns,
                       /*depth=*/0, max_depth, fGeoHolder);
    }

    std::cout << "[EventDisplay] loaded " << fGeoHolder->NumChildren()
              << " geometry shapes\n";

    if (fAutoCenterGeo) AutoCenterGeometry();
}

// -----------------------------------------------------------------------------
void EventDisplay::AutoCenterGeometry()
{
    if (!fGeoHolder || fGeoHolder->NumChildren() == 0) return;

    // Combined bounding box, in scene units.  We use each shape's
    // translation point AND its half-extents, so a single very large
    // shape doesn't fool us into thinking everything is at the origin.
    double xmin = std::numeric_limits<double>::infinity();
    double xmax = -xmin, ymin = xmin, ymax = -xmin, zmin = xmin, zmax = -xmin;

    for (auto* child : fGeoHolder->RefChildren()) {
        auto* gs = dynamic_cast<REX::REveGeoShape*>(child);
        if (!gs) continue;

        // Get translation via REveTrans::GetPos(x, y, z).
        double tx, ty, tz;
        gs->RefMainTrans().GetPos(tx, ty, tz);

        // Half-extents: TGeoBBox is the base of nearly all concrete shape
        // types (TGeoTube, TGeoSphere, TGeoCompositeShape, ...), so a
        // dynamic_cast catches the common cases.  For an exotic shape
        // that doesn't derive from TGeoBBox we just use the translation
        // point alone -- close enough for centering.
        TGeoShape* sh = gs->GetShape();
        TGeoBBox*  bb = dynamic_cast<TGeoBBox*>(sh);
        const double dx = bb ? bb->GetDX() : 0.0;
        const double dy = bb ? bb->GetDY() : 0.0;
        const double dz = bb ? bb->GetDZ() : 0.0;

        xmin = std::min(xmin, tx - dx); xmax = std::max(xmax, tx + dx);
        ymin = std::min(ymin, ty - dy); ymax = std::max(ymax, ty + dy);
        zmin = std::min(zmin, tz - dz); zmax = std::max(zmax, tz + dz);
    }

    fGeoOff[0] = 0.5 * (xmin + xmax);
    fGeoOff[1] = 0.5 * (ymin + ymax);
    fGeoOff[2] = 0.5 * (zmin + zmax);

    // Apply: subtract offset from each shape's translation.
    for (auto* child : fGeoHolder->RefChildren()) {
        auto* gs = dynamic_cast<REX::REveGeoShape*>(child);
        if (!gs) continue;
        REX::REveTrans& tr = gs->RefMainTrans();
        double tx, ty, tz;
        tr.GetPos(tx, ty, tz);
        tr.SetPos(tx - fGeoOff[0], ty - fGeoOff[1], tz - fGeoOff[2]);
    }

    std::cout << "[EventDisplay] auto-centered detector. geo offset (scene units) = ("
              << fGeoOff[0] << ", " << fGeoOff[1] << ", " << fGeoOff[2] << ")\n"
              << "                bbox before shift: x[" << xmin << ", " << xmax
              << "] y[" << ymin << ", " << ymax
              << "] z[" << zmin << ", " << zmax << "]\n"
              << "                bbox after  shift: x["
              << xmin - fGeoOff[0] << ", " << xmax - fGeoOff[0] << "] y["
              << ymin - fGeoOff[1] << ", " << ymax - fGeoOff[1] << "] z["
              << zmin - fGeoOff[2] << ", " << zmax - fGeoOff[2] << "]\n";

    // Remember scene-space bbox so hit auto-centring can align to it.
    fGeoSceneMin[0] = xmin - fGeoOff[0]; fGeoSceneMax[0] = xmax - fGeoOff[0];
    fGeoSceneMin[1] = ymin - fGeoOff[1]; fGeoSceneMax[1] = ymax - fGeoOff[1];
    fGeoSceneMin[2] = zmin - fGeoOff[2]; fGeoSceneMax[2] = zmax - fGeoOff[2];
    fGeoBboxValid = true;

    fAutoCenterGeo = false;   // mark as already done
}

// -----------------------------------------------------------------------------
void EventDisplay::AutoCenterOnFirstEvent()
{
    if (!fTree || fTree->GetEntries() == 0) return;

    // Pull entry 0 -- branch addresses are already wired up by OpenEventFile.
    fTree->GetEntry(0);
    if (!fXg || fXg->empty()) return;

    double xmin = std::numeric_limits<double>::infinity();
    double xmax = -xmin, ymin = xmin, ymax = -xmin, zmin = xmin, zmax = -xmin;
    for (size_t k = 0; k < fXg->size(); ++k) {
        const double x = (*fXg)[k] * fHitScale;
        const double y = (*fYg)[k] * fHitScale;
        const double z = (*fZg)[k] * fHitScale;
        xmin = std::min(xmin, x); xmax = std::max(xmax, x);
        ymin = std::min(ymin, y); ymax = std::max(ymax, y);
        zmin = std::min(zmin, z); zmax = std::max(zmax, z);
    }

    // x, y: center on the hit bbox centre (the shower lands near the
    // beam axis, which is what the user's eye expects to see at the
    // middle of the calorimeter face).
    fHitOff[0] = 0.5 * (xmin + xmax);
    fHitOff[1] = 0.5 * (ymin + ymax);

    // z: if the geometry is loaded and we know its scene bbox, align
    // the hits' front face to the geometry's front face -- this is the
    // physically meaningful choice for a beam-test calorimeter.
    // Otherwise, fall back to centring on the hits' bbox.
    if (fGeoBboxValid && fAlignHitsZFront) {
        // Want: hit_z_min_in_scene == fGeoSceneMin[2]
        //   <=>  zmin - fHitOff[2]   == fGeoSceneMin[2]
        //   <=>  fHitOff[2]          == zmin - fGeoSceneMin[2]
        fHitOff[2] = zmin - fGeoSceneMin[2];
    } else {
        fHitOff[2] = 0.5 * (zmin + zmax);
    }

    fAutoCenterHits = false;   // mark as already done

    std::cout << "[EventDisplay] auto-centered on first event hits. hit offset = ("
              << fHitOff[0] << ", " << fHitOff[1] << ", " << fHitOff[2] << ")\n"
              << "                hit bbox before shift: x[" << xmin << ", " << xmax
              << "] y[" << ymin << ", " << ymax
              << "] z[" << zmin << ", " << zmax << "]\n";
    if (fGeoBboxValid && fAlignHitsZFront) {
        std::cout << "                z aligned: hit front face -> geo front face ("
                  << fGeoSceneMin[2] << ")\n";
    }
}

// -----------------------------------------------------------------------------
void EventDisplay::BuildGeoShapes(TGeoNode* node,
                                  const TGeoHMatrix& parent_mtx,
                                  const std::vector<std::regex>& patterns,
                                  int depth, int max_depth,
                                  REX::REveElement* parent)
{
    if (!node) return;
    if (max_depth >= 0 && depth > max_depth) return;

    // Cumulative global transform = parent * local
    TGeoHMatrix global(parent_mtx);
    if (node->GetMatrix()) global.Multiply(node->GetMatrix());

    TGeoVolume* vol = node->GetVolume();
    if (!vol) return;
    const std::string vname = vol->GetName();

    bool matched = false;
    for (const auto& re : patterns) {
        if (std::regex_search(vname, re)) { matched = true; break; }
    }

    if (matched) {
        TGeoShape* shape = vol->GetShape();
        if (shape) {
            auto* gs = new REX::REveGeoShape(vname.c_str());
            gs->SetShape(shape);
            gs->RefMainTrans().SetFrom(global);

            // Pick a colour AND a transparency by guessing from the
            // volume name.  Absorbers (Lead, Iron) get high transparency
            // and a light tone so a stacked column of them doesn't
            // accumulate into a wall when viewed edge-on; scintillators
            // are the layers you actually want to see, so they're
            // brighter and more solid.
            int col   = kAzure - 9;
            int trans = 80;
            if      (vname.find("Lead")    != std::string::npos) { col = kGray;        trans = 92; }
            else if (vname.find("Iron")    != std::string::npos) { col = kGray + 1;    trans = 92; }
            else if (vname.find("WidePVT") != std::string::npos) { col = kAzure + 1;   trans = 55; }
            else if (vname.find("ThinPS")  != std::string::npos) { col = kOrange + 1;  trans = 55; }
            else if (vname.find("HPL")     != std::string::npos) { col = kGreen + 2;   trans = 75; }

            gs->SetMainColor(col);
            gs->SetMainTransparency(trans);
            gs->SetNSegments(40);

            parent->AddElement(gs);
        }
        // Don't descend into matched subtree.  Comment out the return if you
        // want to also see the children (e.g. individual scintillator strips).
        return;
    }

    int nd = node->GetNdaughters();
    for (int i = 0; i < nd; ++i) {
        BuildGeoShapes(node->GetDaughter(i), global, patterns,
                       depth + 1, max_depth, parent);
    }
}

// -----------------------------------------------------------------------------
void EventDisplay::OpenEventFile(const std::string& root_path,
                                 const std::string& tree_name)
{
    fFile = TFile::Open(root_path.c_str(), "READ");
    if (!fFile || fFile->IsZombie()) {
        std::cerr << "[EventDisplay] cannot open " << root_path << "\n";
        return;
    }
    fTree = static_cast<TTree*>(fFile->Get(tree_name.c_str()));
    if (!fTree) {
        std::cerr << "[EventDisplay] tree '" << tree_name << "' not found\n";
        return;
    }

    fTree->SetBranchAddress("edep",     &fEdep);
    fTree->SetBranchAddress("x_global", &fXg);
    fTree->SetBranchAddress("y_global", &fYg);
    fTree->SetBranchAddress("z_global", &fZg);
    fTree->SetBranchAddress("type",     &fType);
    fTree->SetBranchAddress("section",  &fSection);
    fTree->SetBranchAddress("layer",    &fLayer);
    fTree->SetBranchAddress("hcal",     &fHcal);
    fTree->SetBranchAddress("hexant",   &fHexant);

    std::cout << "[EventDisplay] " << fTree->GetEntries()
              << " events in " << root_path << "\n";

    // Always auto-centre hits on their own bbox (independent of any
    // geometry centring), so they show up at the scene origin even when
    // the simulation's coordinate frame doesn't match the GDML's.
    if (fAutoCenterHits) {
        AutoCenterOnFirstEvent();
    }
}

// -----------------------------------------------------------------------------
void EventDisplay::GotoEvent(Long64_t i)
{
    if (!fTree) return;
    if (i < 0 || i >= fTree->GetEntries()) {
        std::cerr << "[EventDisplay] event " << i << " out of range\n";
        return;
    }

    fTree->GetEntry(i);
    fCurrent = i;

    REX::REveManager::ChangeGuard guard;

    // Clear previous event's elements
    fEventHolder->DestroyElements();

    if (!fEdep || fEdep->empty()) {
        std::cout << "[EventDisplay] event " << i << " has no hits\n";
        return;
    }

    // Energy range (used only for the log-binning of point-set bins).
    const double emin = *std::min_element(fEdep->begin(), fEdep->end());
    const double emax = std::max(emin + 1e-6,
                                 *std::max_element(fEdep->begin(), fEdep->end()));

    // ---------------------------------------------------------------------
    // Hit rendering via REvePointSet bins (one PointSet per energy band).
    // REvePointSet is the simplest, most stable element in REve, so this
    // path bypasses any rendering quirk REveBoxSet might have on the local
    // ROOT build.  Once you can see these points, you can swap back to
    // boxes for a more "calorimetric" look.
    //
    // Energy is split into 5 log bins because edep spans many orders of
    // magnitude (sub-MeV MIP-like deposits up to multi-100-MeV showers).
    // ---------------------------------------------------------------------
    constexpr int kNBins = 5;
    const Color_t bin_colors[kNBins] = {
        kAzure + 1, kCyan + 1, kGreen + 1, kOrange + 7, kRed + 1
    };
    const Float_t bin_sizes[kNBins] = { 4.5f, 5.5f, 6.5f, 7.5f, 8.5f };

    REX::REvePointSet* pbin[kNBins];
    for (int b = 0; b < kNBins; ++b) {
        pbin[b] = new REX::REvePointSet(
            Form("hits_bin_%d", b),
            Form("Energy bin %d", b));
        pbin[b]->SetMarkerStyle(20);            // filled circle
        pbin[b]->SetMarkerSize(bin_sizes[b]);
        pbin[b]->SetMarkerColor(bin_colors[b]);
    }

    // Log-binning of edep.  Guard for emin == 0 (zero-edep hits get pushed
    // to bin 0 by the clamp below).
    const double log_lo = std::log10(std::max(emin, 1e-9));
    const double log_hi = std::log10(std::max(emax, log_lo * 10));
    const double log_w  = std::max(log_hi - log_lo, 1e-6);

    // Spatial bbox of hits, in *scene* units (after hit-scale and offset).
    double xmin = std::numeric_limits<double>::infinity(), xmax = -xmin;
    double ymin = xmin, ymax = -xmin;
    double zmin = xmin, zmax = -xmin;

    for (size_t k = 0; k < fEdep->size(); ++k) {
        const float x = static_cast<float>((*fXg)[k] * fHitScale - fHitOff[0]);
        const float y = static_cast<float>((*fYg)[k] * fHitScale - fHitOff[1]);
        const float z = static_cast<float>((*fZg)[k] * fHitScale - fHitOff[2]);

        const double le = std::log10(std::max((*fEdep)[k], 1e-9));
        int b = static_cast<int>(((le - log_lo) / log_w) * kNBins);
        if (b < 0) b = 0;
        if (b >= kNBins) b = kNBins - 1;
        pbin[b]->SetNextPoint(x, y, z);

        xmin = std::min<double>(xmin, x); xmax = std::max<double>(xmax, x);
        ymin = std::min<double>(ymin, y); ymax = std::max<double>(ymax, y);
        zmin = std::min<double>(zmin, z); zmax = std::max<double>(zmax, z);
    }

    // Ensure at least one always-visible reference point at the scene
    // origin -- handy for verifying that the viewer is rendering at all.
    auto* anchor = new REX::REvePointSet("anchor", "scene origin");
    anchor->SetMarkerStyle(2);   // big "+"
    anchor->SetMarkerSize(3.0);
    anchor->SetMarkerColor(kWhite);
    anchor->SetNextPoint(0.0, 0.0, 0.0);
    fEventHolder->AddElement(anchor);

    for (int b = 0; b < kNBins; ++b) {
        fEventHolder->AddElement(pbin[b]);
    }

    // ---- Optional: a per-event "cluster" per (section, layer) bucket ------
    // The user mentioned "something akin to clusters".  As a placeholder we
    // simply group hits by (section, layer) and report multiplicity / total
    // edep in element titles -- you'll want to replace this with a real
    // clustering algorithm (nearest-neighbour, topo-clustering, etc.).
    //
    // To turn a group into a *visible* cluster element, add a second
    // REveBoxSet with the cluster centroids, or a REvePointSet.
    // ----------------------------------------------------------------------

    std::cout << "[EventDisplay] event " << i << ": "
              << fEdep->size() << " hits, edep \u2208 ["
              << emin << ", " << emax << "] MeV\n"
              << "                hit bbox (after scale + offset): x[" << xmin << ", " << xmax
              << "] y[" << ymin << ", " << ymax
              << "] z[" << zmin << ", " << zmax << "]\n";
}

// -----------------------------------------------------------------------------
//                              main()
// -----------------------------------------------------------------------------
namespace {
    EventDisplay* g_display = nullptr;
}
// These free functions are exposed to the REve C++-interpreter UI so you
// can wire buttons to them. From the browser console / command line you
// can also type:   eve_next();   eve_prev();   eve_goto(42);
extern "C" void eve_next() { if (g_display) g_display->NextEvent(); }
extern "C" void eve_prev() { if (g_display) g_display->PrevEvent(); }
extern "C" void eve_goto(Long64_t i) { if (g_display) g_display->GotoEvent(i); }

int main(int argc, char** argv)
{
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <detector.gdml> <events.root> [tree_name]\n";
        return 1;
    }

    // Save argv contents BEFORE constructing TApplication.  Its constructor
    // takes &argc, &argv and (in some ROOT versions) clobbers entries it
    // doesn't recognise -- which is why the geometry path showed up empty.
    const std::string gdml_path = argv[1];
    const std::string root_path = argv[2];
    const std::string tree_name = (argc >= 4) ? argv[3] : "calo_events";

    // ROOT needs a TApplication to drive the web-server event loop.
    TApplication app("eve_app", &argc, argv);

    EventDisplay ed;
    g_display = &ed;
    ed.Init();

    // ---- Geometry filters -------------------------------------------------
    // Keep:
    //   - the 86 ECAL layer envelopes (Lead + WidePVT + ThinPS) across
    //     all modules (MX1Y1, MX1Y2, MX2Y1, MX2Y2, MX2Y3, ...);
    //   - the HCAL scintillator layer envelopes (HCAL_GLn_...);
    //   - the iron absorber plates (IronPlateLog).
    // Skip the thousands of HPL_* hodoscope fibres -- they would tank
    // performance.  Edit this list to taste.
    const std::vector<std::string> include = {
        R"(^ECAL_GL\d+_)",                // ECAL Lead + scintillator layers
        R"(^HCAL_GL\d+_)",                // HCAL scintillator layers
        R"(^IronPlate)",                  // HCAL iron absorbers
        // R"(^HPL_FiberCoreLog)",        // uncomment to also draw fibre cores
    };
    ed.LoadGeometry(gdml_path, include);

    ed.OpenEventFile(root_path, tree_name);

    // Geometry and hit data may live in different absolute world frames
    // in this dataset (geometry centred near z=37 m, hits near z=96 m).
    // Since we keep separate offsets for the two, both end up at the
    // scene origin and are visually overlapping in the viewer -- handy
    // for inspection, but be aware: their physical positions in the
    // world differ by ~60 m, so this overlay is a *visualisation*, not
    // physical ground truth.

    if (ed.NumEvents() > 0) ed.GotoEvent(0);

    ed.Manager()->Show();   // launches default browser at the local URL

    app.Run();
    return 0;
}
