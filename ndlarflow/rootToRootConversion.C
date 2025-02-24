// ROOT includes
#include "TROOT.h"
#include "TInterpreter.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

// C++ includes
#include <iostream>
#include <vector>

// NB: this is not so elegant for a first pass, but it seems to work

void rootToRootConversion(
  const bool isMC=true,
  const std::string fname="MicroProdN3p1_NDLAr_2E18_FHC.flow.nu.0000001.FLOW.hdf5_hits_uproot.root",
  const std::string outname="MicroProdN3p1_NDLAr_2E18_FHC.flow.nu.0000001.FLOW.hdf5_hits.root" )
{
    // generate the dictionaries for vectors of vectors : see e.g.
    //   https://root-forum.cern.ch/t/problem-in-accessing-vector-vector-float/27983
    // though others also do this.
    gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
    gInterpreter->GenerateDictionary("vector<vector<long> >", "vector");
    gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");

    TFile *f = new TFile(fname.c_str(),"read");
    const int MaxDepthArray = 10000;
    const int MaxDepthArrayMC = isMC ? MaxDepthArray : 1;
    const int MaxDepthArrayNu = isMC ? 500 : 1;

    // Get the data products FROM the tree
    int in_run, in_subrun, in_event, in_subevent, in_event_start_t, in_event_end_t, in_unix_ts;

    // hits
    int Nhits;
    unsigned short matches[MaxDepthArray];
    Float_t in_x[MaxDepthArray];
    Float_t in_y[MaxDepthArray];
    Float_t in_z[MaxDepthArray];
    Float_t in_ts[MaxDepthArray];
    Float_t in_charge[MaxDepthArray];
    Float_t in_E[MaxDepthArray];

    // Hit truth (currently cappped at 10 per hit)
    int     nmatchesinsubevent(0);
    Float_t packetFrac[MaxDepthArrayMC];
    Long_t  particleID[MaxDepthArrayMC];
    Long_t  particleIDLocal[MaxDepthArrayMC];
    Int_t   pdg[MaxDepthArrayMC];
    Long_t  vertexID[MaxDepthArrayMC];
    Long_t  segmentID[MaxDepthArrayMC];

    // MCP
    // TODO: we should probably save the max number of MCParticles and neutrinos as well, for now we'll assume 10*NMaxHits which is hopefully okay but not guaranteed
    int     nmcpinsubevent(0);
    Float_t in_mcp_energy[MaxDepthArrayMC];
    Float_t in_mcp_px[MaxDepthArrayMC];
    Float_t in_mcp_py[MaxDepthArrayMC];
    Float_t in_mcp_pz[MaxDepthArrayMC];
    Float_t in_mcp_startx[MaxDepthArrayMC];
    Float_t in_mcp_starty[MaxDepthArrayMC];
    Float_t in_mcp_startz[MaxDepthArrayMC];
    Float_t in_mcp_endx[MaxDepthArrayMC];
    Float_t in_mcp_endy[MaxDepthArrayMC];
    Float_t in_mcp_endz[MaxDepthArrayMC];
    Int_t   in_mcp_pdg[MaxDepthArrayMC];
    Long_t  in_mcp_nuid[MaxDepthArrayMC];
    Long_t  in_mcp_vertex_id[MaxDepthArrayMC];
    Long_t  in_mcp_idLocal[MaxDepthArrayMC];
    Long_t  in_mcp_id[MaxDepthArrayMC];
    Long_t  in_mcp_mother[MaxDepthArrayMC];
    // Neutrinos: see caveat above, but here we'll have a hard cap at 500
    int     nnuinsubevent(0);
    Long_t  in_nuID[MaxDepthArrayNu];
    Int_t   in_nuPDG[MaxDepthArrayNu];
    Int_t   in_mode[MaxDepthArrayNu];
    Int_t   in_ccnc[MaxDepthArrayNu];
    Float_t in_nue[MaxDepthArrayNu];
    Float_t in_nupx[MaxDepthArrayNu];
    Float_t in_nupy[MaxDepthArrayNu];
    Float_t in_nupz[MaxDepthArrayNu];
    Float_t in_nuvtxx[MaxDepthArrayNu];
    Float_t in_nuvtxy[MaxDepthArrayNu];
    Float_t in_nuvtxz[MaxDepthArrayNu];

    TTree* tr = (TTree*)f->Get("events");
    tr->SetBranchAddress("run",&in_run);
    tr->SetBranchAddress("subrun",&in_subrun);
    tr->SetBranchAddress("event",&in_event);
    tr->SetBranchAddress("unix_ts",&in_unix_ts);
    tr->SetBranchAddress("event_start_t",&in_event_start_t);
    tr->SetBranchAddress("event_end_t",&in_event_end_t);
    tr->SetBranchAddress("subevent",&in_subevent);

    tr->SetBranchAddress("nx",&Nhits);
    tr->SetBranchAddress("x",&in_x);
    tr->SetBranchAddress("y",&in_y);
    tr->SetBranchAddress("z",&in_z);
    tr->SetBranchAddress("ts",&in_ts);
    tr->SetBranchAddress("charge",&in_charge);
    tr->SetBranchAddress("E",&in_E);

    if ( isMC ) {
        tr->SetBranchAddress("matches",&matches);
        tr->SetBranchAddress("hit_packetFrac",&packetFrac);
        tr->SetBranchAddress("hit_particleID",&particleID);
        tr->SetBranchAddress("hit_particleIDLocal",&particleIDLocal);
        tr->SetBranchAddress("hit_pdg",&pdg);
        tr->SetBranchAddress("hit_vertexID",&vertexID);
        tr->SetBranchAddress("hit_segmentID",&segmentID);
        tr->SetBranchAddress("nhit_packetFrac",&nmatchesinsubevent);

        tr->SetBranchAddress("mcp_energy", &in_mcp_energy);
        tr->SetBranchAddress("mcp_px", &in_mcp_px);
        tr->SetBranchAddress("mcp_py", &in_mcp_py);
        tr->SetBranchAddress("mcp_pz", &in_mcp_pz);
        tr->SetBranchAddress("mcp_startx", &in_mcp_startx);
        tr->SetBranchAddress("mcp_starty", &in_mcp_starty);
        tr->SetBranchAddress("mcp_startz", &in_mcp_startz);
        tr->SetBranchAddress("mcp_endx", &in_mcp_endx);
        tr->SetBranchAddress("mcp_endy", &in_mcp_endy);
        tr->SetBranchAddress("mcp_endz", &in_mcp_endz);
        tr->SetBranchAddress("mcp_pdg", &in_mcp_pdg);
        tr->SetBranchAddress("mcp_nuid", &in_mcp_nuid);
        tr->SetBranchAddress("mcp_vertex_id", &in_mcp_vertex_id);
        tr->SetBranchAddress("mcp_idLocal", &in_mcp_idLocal);
        tr->SetBranchAddress("mcp_id", &in_mcp_id);
        tr->SetBranchAddress("mcp_mother", &in_mcp_mother);
        tr->SetBranchAddress("nmcp_energy",&nmcpinsubevent);

        tr->SetBranchAddress("nuID", &in_nuID); // both nuID and vertex_id appear to use this value...
        tr->SetBranchAddress("nuPDG", &in_nuPDG);
        tr->SetBranchAddress("mode", &in_mode);
        tr->SetBranchAddress("ccnc", &in_ccnc);
        tr->SetBranchAddress("nue", &in_nue);
        tr->SetBranchAddress("nupx", &in_nupx);
        tr->SetBranchAddress("nupy", &in_nupy);
        tr->SetBranchAddress("nupz", &in_nupz);
        tr->SetBranchAddress("nuvtxx", &in_nuvtxx);
        tr->SetBranchAddress("nuvtxy", &in_nuvtxy);
        tr->SetBranchAddress("nuvtxz", &in_nuvtxz);
        tr->SetBranchAddress("nnuID",&nnuinsubevent);
    }

    long sum=0;

    long NEvents = tr->GetEntries();

    std::cout << "Loaded in tree with " << NEvents << " entries." << std::endl;

    // OUTPUT TREE
    int run, subrun, event, event_start_t, event_end_t, unix_ts;

    std::vector<float> x;
    std::vector<float> y;
    std::vector<float> z;
    std::vector<float> ts;
    std::vector<float> E;
    std::vector<float> charge;

    // segmentIndex, particleIndex unfilled!
    std::vector< std::vector<int> >   hit_pdg;
    std::vector< std::vector<long> >  hit_segmentID;
    std::vector< std::vector<long> >  hit_particleID;
    std::vector< std::vector<long> >  hit_particleIDLocal;
    std::vector< std::vector<long> >  hit_vertexID;
    std::vector< std::vector<float> > hit_packetFrac;
    
    std::vector<float> mcp_px;
    std::vector<float> mcp_py;
    std::vector<float> mcp_pz;
    std::vector<long> mcp_id;
    std::vector<long> mcp_idLocal;
    std::vector<long> mcp_nuid;
    std::vector<long> mcp_vertex_id;
    std::vector<int> mcp_pdg;
    std::vector<long> mcp_mother;
    std::vector<float> mcp_energy;
    std::vector<float> mcp_startx;
    std::vector<float> mcp_starty;
    std::vector<float> mcp_startz;
    std::vector<float> mcp_endx;
    std::vector<float> mcp_endy;
    std::vector<float> mcp_endz;

    std::vector<float> nuvtxx;
    std::vector<float> nuvtxy;
    std::vector<float> nuvtxz;
    std::vector<float> nupx;
    std::vector<float> nupy;
    std::vector<float> nupz;
    std::vector<float> nue;
    std::vector<long> nuID;
    std::vector<long> vertex_id;
    std::vector<int> nuPDG;
    std::vector<int> mode;
    std::vector<int> ccnc;

    TFile *outFile = new TFile(outname.c_str(),"recreate");

    TTree *outgoingTree = new TTree("events", "events");
    outgoingTree->Branch("run", &run);
    outgoingTree->Branch("subrun", &subrun);
    outgoingTree->Branch("event", &event);
    outgoingTree->Branch("event_start_t", &event_start_t);
    outgoingTree->Branch("event_end_t", &event_end_t);
    outgoingTree->Branch("unix_ts", &unix_ts);
    outgoingTree->Branch("x", &x);
    outgoingTree->Branch("y", &y);
    outgoingTree->Branch("z", &z);
    outgoingTree->Branch("ts", &ts);
    outgoingTree->Branch("E", &E);
    outgoingTree->Branch("charge", &charge);
    if ( isMC ) {
        outgoingTree->Branch("hit_pdg", &hit_pdg);
        outgoingTree->Branch("hit_segmentID", &hit_segmentID);
        outgoingTree->Branch("hit_particleID", &hit_particleID);
        outgoingTree->Branch("hit_particleIDLocal", &hit_particleIDLocal);
        outgoingTree->Branch("hit_vertexID", &hit_vertexID);
        outgoingTree->Branch("hit_packetFrac", &hit_packetFrac);
        outgoingTree->Branch("mcp_px", &mcp_px);
        outgoingTree->Branch("mcp_py", &mcp_py);
        outgoingTree->Branch("mcp_pz", &mcp_pz);
        outgoingTree->Branch("mcp_id", &mcp_id);
        outgoingTree->Branch("mcp_idLocal", &mcp_idLocal);
        outgoingTree->Branch("mcp_nuid", &mcp_nuid);
        outgoingTree->Branch("mcp_vertex_id", &mcp_vertex_id);
        outgoingTree->Branch("mcp_pdg", &mcp_pdg);
        outgoingTree->Branch("mcp_mother", &mcp_mother);
        outgoingTree->Branch("mcp_energy", &mcp_energy);
        outgoingTree->Branch("mcp_startx", &mcp_startx);
        outgoingTree->Branch("mcp_starty", &mcp_starty);
        outgoingTree->Branch("mcp_startz", &mcp_startz);
        outgoingTree->Branch("mcp_endx", &mcp_endx);
        outgoingTree->Branch("mcp_endy", &mcp_endy);
        outgoingTree->Branch("mcp_endz", &mcp_endz);
        outgoingTree->Branch("nuvtxx", &nuvtxx);
        outgoingTree->Branch("nuvtxy", &nuvtxy);
        outgoingTree->Branch("nuvtxz", &nuvtxz);
        outgoingTree->Branch("nupx", &nupx);
        outgoingTree->Branch("nupy", &nupy);
        outgoingTree->Branch("nupz", &nupz);
        outgoingTree->Branch("nue", &nue);
        outgoingTree->Branch("nuID", &nuID);
        outgoingTree->Branch("vertex_id", &vertex_id);
        outgoingTree->Branch("nuPDG", &nuPDG);
        outgoingTree->Branch("mode", &mode);
        outgoingTree->Branch("ccnc", &ccnc);
    }

    int thisRun=0;
    int thisSubRun=0;
    int thisEvent=0;

    unsigned long sum_matches=0;
    std::vector<short> all_matches;
    std::vector<int>   all_hit_pdg;
    std::vector<long>  all_hit_segmentID;
    std::vector<long>  all_hit_particleID;
    std::vector<long>  all_hit_particleIDLocal;
    std::vector<long>  all_hit_vertexID;
    std::vector<float> all_hit_packetFrac;

    for (unsigned int idx=0; idx<=NEvents; ++idx) {
        std::cout << in_run << ":" << in_subrun << ":" << in_event << ":" << in_subevent << std::endl;
        if (idx!=NEvents){
            tr->GetEntry(idx);
        }
        if ( idx==0 ){
            thisRun = in_run;
            thisSubRun = in_subrun;
            thisEvent = in_event;
        }
        if ( idx > 0 && ( (in_run!=thisRun || in_subrun!=thisSubRun || in_event!=thisEvent) || idx==NEvents ) ) {
            if ( isMC ) {
                // Something has changed... finish with the matches, fill the tree, and reset
                // Vectors of matches
                if ( sum_matches!=all_hit_packetFrac.size() ){
                    std::cout << "WARNING!!! For R" << thisRun << ", S" << thisSubRun << ", E"
                            << thisEvent << " -- the summed number of matches expected (" << sum_matches << ") "
                            << "does **NOT** match the size of the packetFrac vector (" << all_hit_packetFrac.size()
                            << "). This is a sign of trouble." << std::endl;
                }
                unsigned long matchIndex = 0;
                std::vector<int>   this_hit_pdg;
                std::vector<long>  this_hit_segmentID;
                std::vector<long>  this_hit_particleID;
                std::vector<long>  this_hit_particleIDLocal;
                std::vector<long>  this_hit_vertexID;
                std::vector<float> this_hit_packetFrac;
                for ( unsigned int idxMatch=0; idxMatch<=(unsigned int)all_hit_packetFrac.size(); ++idxMatch ){
                    if( this_hit_pdg.size() == (unsigned int)all_matches[matchIndex] || idxMatch==(unsigned int)all_hit_packetFrac.size()) {
                        hit_pdg.push_back( this_hit_pdg );
                        hit_segmentID.push_back( this_hit_segmentID );
                        hit_particleID.push_back( this_hit_particleID );
                        hit_particleIDLocal.push_back( this_hit_particleIDLocal );
                        hit_vertexID.push_back( this_hit_vertexID );
                        hit_packetFrac.push_back( this_hit_packetFrac );
                        // mini-reset
                        matchIndex+=1;
                        this_hit_pdg.clear();
                        this_hit_segmentID.clear();
                        this_hit_particleID.clear();
                        this_hit_particleIDLocal.clear();
                        this_hit_vertexID.clear();
                        this_hit_packetFrac.clear();
                        // break if == end
                        if ( idxMatch == all_hit_packetFrac.size() ) break;
                    }
                    this_hit_pdg.push_back(all_hit_pdg[idxMatch]);
                    this_hit_segmentID.push_back(all_hit_segmentID[idxMatch]);
                    this_hit_particleID.push_back(all_hit_particleID[idxMatch]);
                    this_hit_particleIDLocal.push_back(all_hit_particleIDLocal[idxMatch]);
                    this_hit_vertexID.push_back(all_hit_vertexID[idxMatch]);
                    this_hit_packetFrac.push_back(all_hit_packetFrac[idxMatch]);
                }
            }
            // Fill
            outgoingTree->Fill();
            // Clear
            x.clear();
            y.clear();
            z.clear();
            ts.clear();
            E.clear();
            charge.clear();
            hit_pdg.clear();
            hit_segmentID.clear();
            hit_particleID.clear();
            hit_particleIDLocal.clear();
            hit_vertexID.clear();
            hit_packetFrac.clear();
            mcp_px.clear();
            mcp_py.clear();
            mcp_pz.clear();
            mcp_id.clear();
            mcp_idLocal.clear();
            mcp_nuid.clear();
            mcp_vertex_id.clear();
            mcp_pdg.clear();
            mcp_mother.clear();
            mcp_energy.clear();
            mcp_startx.clear();
            mcp_starty.clear();
            mcp_startz.clear();
            mcp_endx.clear();
            mcp_endy.clear();
            mcp_endz.clear();
            nuvtxx.clear();
            nuvtxy.clear();
            nuvtxz.clear();
            nupx.clear();
            nupy.clear();
            nupz.clear();
            nue.clear();
            nuID.clear();
            vertex_id.clear();
            nuPDG.clear();
            mode.clear();
            ccnc.clear();
            
            sum_matches=0;
            all_matches.clear();
            all_hit_pdg.clear();
            all_hit_segmentID.clear();
            all_hit_particleID.clear();
            all_hit_particleIDLocal.clear();
            all_hit_vertexID.clear();
            all_hit_packetFrac.clear();

            // Set up the current run/subrun/event to look for
            thisRun = in_run;
            thisSubRun = in_subrun;
            thisEvent = in_event;
            
            if (idx==NEvents) break; // We're done!
        }

        run = in_run;
        subrun = in_subrun;
        event = in_event;
        event_start_t = in_event_start_t;
        event_end_t = in_event_end_t;
        unix_ts = in_unix_ts;
        
        // fill up the vectors for as much stuff as we can in this subevent:
        for ( unsigned int idxHit=0; idxHit<(unsigned int)Nhits; ++idxHit ) {
            x.push_back(in_x[idxHit]);
            y.push_back(in_y[idxHit]);
            z.push_back(in_z[idxHit]);
            ts.push_back(in_ts[idxHit]);
            E.push_back(in_E[idxHit]);
            charge.push_back(in_charge[idxHit]);
            if ( isMC ) {
                all_matches.push_back(matches[idxHit]);
                sum_matches+=(unsigned long)matches[idxHit];
            }
        }

        if ( isMC ) {
            // Because the arrays of matches and hits need not be sync'ed, we'll fill up
            //   vectors will all the info for the event here and then break it into sub-
            //   vectors at the end of the event...
            for ( unsigned int idxMatch=0; idxMatch<(unsigned int)nmatchesinsubevent; ++idxMatch ) {
                all_hit_pdg.push_back(pdg[idxMatch]);
                all_hit_segmentID.push_back(segmentID[idxMatch]);
                all_hit_particleID.push_back(particleID[idxMatch]);
                all_hit_particleIDLocal.push_back(particleIDLocal[idxMatch]);
                all_hit_vertexID.push_back(vertexID[idxMatch]);
                all_hit_packetFrac.push_back(packetFrac[idxMatch]);
            }

            // MCParticles
            for ( unsigned int idxMCPart=0; idxMCPart<(unsigned int)nmcpinsubevent; ++idxMCPart ) {
                mcp_px.push_back(in_mcp_px[idxMCPart]);
                mcp_py.push_back(in_mcp_py[idxMCPart]);
                mcp_pz.push_back(in_mcp_pz[idxMCPart]);
                mcp_id.push_back(in_mcp_id[idxMCPart]);
                mcp_idLocal.push_back(in_mcp_idLocal[idxMCPart]);
                mcp_nuid.push_back(in_mcp_nuid[idxMCPart]);
                mcp_vertex_id.push_back(in_mcp_vertex_id[idxMCPart]);
                mcp_pdg.push_back(in_mcp_pdg[idxMCPart]);
                mcp_mother.push_back(in_mcp_mother[idxMCPart]);
                mcp_energy.push_back(in_mcp_energy[idxMCPart]);
                mcp_startx.push_back(in_mcp_startx[idxMCPart]);
                mcp_starty.push_back(in_mcp_starty[idxMCPart]);
                mcp_startz.push_back(in_mcp_startz[idxMCPart]);
                mcp_endx.push_back(in_mcp_endx[idxMCPart]);
                mcp_endy.push_back(in_mcp_endy[idxMCPart]);
                mcp_endz.push_back(in_mcp_endz[idxMCPart]);
            }

            // Neutrinos
            for ( unsigned int idxNu=0; idxNu<(unsigned int)nnuinsubevent; ++idxNu ) {
                nuvtxx.push_back(in_nuvtxx[idxNu]);
                nuvtxy.push_back(in_nuvtxy[idxNu]);
                nuvtxz.push_back(in_nuvtxz[idxNu]);
                nupx.push_back(in_nupx[idxNu]);
                nupy.push_back(in_nupy[idxNu]);
                nupz.push_back(in_nupz[idxNu]);
                nue.push_back(in_nue[idxNu]);
                nuID.push_back(in_nuID[idxNu]);
                vertex_id.push_back(in_nuID[idxNu]); // appears to be the same as nuID in current h5->ROOT script
                nuPDG.push_back(in_nuPDG[idxNu]);
                mode.push_back(in_mode[idxNu]);
                ccnc.push_back(in_ccnc[idxNu]);
            }
        }
    }

    outgoingTree->Write();
    outFile->Close();

    std::cout << "DONE!" << std::endl;
}
