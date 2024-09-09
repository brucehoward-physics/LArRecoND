# Adapted for 2x2 conversion for Pandora by Richard Diurba
# Code based on Bern Module data root file conversion by Salvatore Davide Porzio and Saba Parsa
# python h5_to_root_ndlarflow.py filename data
# filename is an argument of the filename (does not need quotations)
# data is to toggle the MC information for simulation (0) and data (1), default is simulation
# Requires standard ROOT, h5py, and numpy 
# Also requires h5flow, a package by Peter Madigan (https://github.com/peter-madigan/h5flow)
#To install, git clone into a directory and then run pip install .
#Please use a virtual environment on Fermilab machines (venv)
######
from h5flow.data import dereference
import h5flow
from array import array
from collections import defaultdict
import h5py
import numpy as np
import os
import ROOT
from ROOT import  TFile
import sys
trueXOffset=0 # Offsets if geometry changes
trueYOffset=0#42+268
trueZOffset=0#-1300
def main(argv=None):
    
    # set input files to be loaded
    datapath=""
    files = [
"/pnfs/dune/persistent/users/noeroy/prod/MiniRun5_1E19_RHC/MiniRun5_1E19_RHC.flow.beta2a/FLOW/0000000/MiniRun5_1E19_RHC.flow.0000055.FLOW.hdf5"]


    if (len(sys.argv)>1):
        if (sys.argv[1]!=None):
            files=[str(sys.argv[1])]
    useData=False
    if (len(sys.argv)>2):
        if (int(sys.argv[2])==1):
            useData=True
    output="/exp/dune/data/users/rdiurba/flowToROOT"
    if (len(sys.argv)>3):
        if (str(sys.argv[3])!=None):
            output=str(sys.argv[3])
    fnum=0;
    # load files in list
    for fname in files:

        fnum=fnum+1;
        f = h5py.File(os.path.join(datapath,fname),'r')
        print(' Processing file', fname)
        flow_out=h5flow.data.H5FlowDataManager(os.path.join(datapath,fname),"r")
        #load events
        events    = f['charge/events/data']

        # fill tree
        outfile=fname.rsplit("/",1)[-1]
        outFileName = output+"/"+outfile+"_hits.root"


        #print('output file : ', '' + outFileName )

        output_file = ROOT.TFile((outFileName),"RECREATE")
        output_tree = ROOT.TTree("events", "events")

        
        #set up variables
        # Event informations
        eventID              = array('i',[0])           # event ID [-]
        event_start_t        = array('i',[0])           # event timestamp start [UNITS?]
        event_end_t          = array('i',[0])           # event timestamp end [UNITS?]       
        event_unix_ts = array("l",[0])
        subrun=array("i",[0])
        run= array("i",[0])
        
        # vectors and variables needed for Pandora
        x=ROOT.std.vector('float')();
        y=ROOT.std.vector('float')();
        z=ROOT.std.vector('float')();
        ts=ROOT.std.vector('float')();
        E=ROOT.std.vector("float")();
        charge=ROOT.std.vector('float')();
        hit_pdg=ROOT.std.vector("std::vector<int>")();
        hit_particleID=ROOT.std.vector("std::vector<int>")();
        hit_segmentID=ROOT.std.vector("std::vector<int>")();
        hit_packetFrac=ROOT.std.vector("std::vector<float>")();
        hit_segmentIndex=ROOT.std.vector("std::vector<int>")();
        hit_particleIndex=ROOT.std.vector("std::vector<int>")();
        hit_interactionIndex=ROOT.std.vector("std::vector<int>")();
        x_uncalib=ROOT.std.vector('float')();
        y_uncalib=ROOT.std.vector('float')();
        z_uncalib=ROOT.std.vector('float')();
        ts_uncalib=ROOT.std.vector('float')();
        E_uncalib=ROOT.std.vector("float")();
        charge_uncalib=ROOT.std.vector('float')();
        hit_pdg_uncalib=ROOT.std.vector("std::vector<int>")();
        hit_particleID_uncalib=ROOT.std.vector("std::vector<int>")();
        hit_segmentID_uncalib=ROOT.std.vector("std::vector<int>")();
        hit_packetFrac_uncalib=ROOT.std.vector("std::vector<float>")();
        hit_segmentIndex_uncalib=ROOT.std.vector("std::vector<int>")();
        hit_particleIndex_uncalib=ROOT.std.vector("std::vector<int>")();
        hit_interactionIndex_uncalib=ROOT.std.vector("std::vector<int>")();
        mcp_px=ROOT.std.vector("float")();
        mcp_py=ROOT.std.vector("float")();
        mcp_pz=ROOT.std.vector("float")();
        mcp_id=ROOT.std.vector("int")();
        mcp_nuid=ROOT.std.vector("int")();
        mcp_vertexID=ROOT.std.vector("int")();
        mcp_pdg=ROOT.std.vector("int")();
        mcp_mother=ROOT.std.vector("int")();
        mcp_energy=ROOT.std.vector("float")();
        mcp_startx=ROOT.std.vector("float")();
        mcp_starty=ROOT.std.vector("float")();
        mcp_startz=ROOT.std.vector("float")();
        mcp_endx=ROOT.std.vector("float")();
        mcp_endy=ROOT.std.vector("float")();
        mcp_endz=ROOT.std.vector("float")();
        nuvtxx=ROOT.std.vector("float")();
        nuvtxy=ROOT.std.vector("float")();
        nuvtxz=ROOT.std.vector("float")();
        nupx=ROOT.std.vector("float")();
        nupy=ROOT.std.vector("float")();
        nupz=ROOT.std.vector("float")();
        nue=ROOT.std.vector("float")();
        nuID=ROOT.std.vector("int")();
        nuPDG=ROOT.std.vector("int")();
        mode=ROOT.std.vector("int")();
        ccnc=ROOT.std.vector("int")(); 
        
        # stetup tree for output       
        output_tree.Branch("event"           ,eventID           ,"eventID/I")
        output_tree.Branch("subrun"           ,subrun           ,"subrun/I")
        output_tree.Branch("unix_ts",event_unix_ts,"unix_ts/L");
        output_tree.Branch("run"           ,run           ,"run/I")
        output_tree.Branch("event_start_t"     ,event_start_t     ,"event_start_t/I")     # 32 bit timestamp (2^32-1 = 2.147483647e9)
        output_tree.Branch("event_end_t"       ,event_end_t       ,"event_end_t/I")       # 32 bit timestamp (2^32-1 = 2.147483647e9)

        output_tree.Branch("x",x)
        output_tree.Branch("y",y)
        output_tree.Branch("z",z)
        output_tree.Branch("ts",ts)
        output_tree.Branch("charge"      ,charge)
        output_tree.Branch("E",E)
        output_tree.Branch("hit_segmentID",hit_segmentID)
        output_tree.Branch("hit_segmentIndex",hit_segmentIndex)
        output_tree.Branch("hit_particleID",hit_particleID)
        output_tree.Branch("hit_particleIndex",hit_particleIndex)
        output_tree.Branch("hit_interactionIndex",hit_interactionIndex);
        output_tree.Branch("hit_pdg",hit_pdg)
        output_tree.Branch("hit_packetFrac",hit_packetFrac)

        output_tree.Branch("x_uncalib",x_uncalib)
        output_tree.Branch("y_uncalib",y_uncalib)
        output_tree.Branch("z_uncalib",z_uncalib)
        output_tree.Branch("ts_uncalib",ts_uncalib)
        output_tree.Branch("charge_uncalib"      ,charge_uncalib)
        output_tree.Branch("E_uncalib",E_uncalib)
        output_tree.Branch("hit_segmentID_uncalib",hit_segmentID_uncalib)
        output_tree.Branch("hit_segmentIndex_uncalib",hit_segmentIndex_uncalib)
        output_tree.Branch("hit_particleID_uncalib",hit_particleID_uncalib)
        output_tree.Branch("hit_particleIndex_uncalib",hit_particleIndex_uncalib)
        output_tree.Branch("hit_interactionIndex_uncalib",hit_interactionIndex_uncalib);
        output_tree.Branch("hit_pdg_uncalib",hit_pdg_uncalib)
        output_tree.Branch("hit_packetFrac_uncalib",hit_packetFrac_uncalib)

        output_tree.Branch("mcp_energy",mcp_energy)
        output_tree.Branch("mcp_pdg",mcp_pdg)
        output_tree.Branch("mcp_nuid",mcp_nuid)
       
        output_tree.Branch("mcp_id",mcp_id)
        output_tree.Branch("mcp_px",mcp_px)
        output_tree.Branch("mcp_py",mcp_py)
        output_tree.Branch("mcp_pz",mcp_pz)
        output_tree.Branch("mcp_mother",mcp_mother)
        output_tree.Branch("mcp_startx",mcp_startx)
        output_tree.Branch("mcp_starty",mcp_starty)
        output_tree.Branch("mcp_startz",mcp_startz)
        output_tree.Branch("mcp_endx",mcp_endx)
        output_tree.Branch("mcp_endy",mcp_endy)
        output_tree.Branch("mcp_endz",mcp_endz) 
        output_tree.Branch("nuID",nuID)
        output_tree.Branch("nue",nue)
        output_tree.Branch("nuPDG",nuPDG)
        output_tree.Branch("nupx",nupx)
        output_tree.Branch("nupy",nupy)
        output_tree.Branch("nupz",nupz)
        output_tree.Branch("nuvtxx",nuvtxx)
        output_tree.Branch("nuvtxy",nuvtxy)
        output_tree.Branch("nuvtxz",nuvtxz)
        output_tree.Branch("mode",mode)
        output_tree.Branch("ccnc",ccnc)



        run[0]=int(0)
        subrun[0]=int(0)
        # loop over events
        for ev_index in range(len(events)):
            # refresh variables
            x.clear()
            y.clear()
            z.clear()
            charge.clear()
            E.clear()
            ts.clear()
            hit_segmentID.clear()
            hit_pdg.clear()
            hit_particleID.clear()
            hit_interactionIndex.clear()
            hit_particleIndex.clear()
            hit_packetFrac.clear()
            hit_segmentIndex.clear()

            x_uncalib.clear()
            y_uncalib.clear()
            z_uncalib.clear()
            charge_uncalib.clear()
            E_uncalib.clear()
            ts_uncalib.clear()
            hit_segmentID_uncalib.clear()
            hit_pdg_uncalib.clear()
            hit_particleID_uncalib.clear()
            hit_interactionIndex_uncalib.clear()
            hit_particleIndex_uncalib.clear()
            hit_packetFrac_uncalib.clear()
            hit_segmentIndex_uncalib.clear()

            mcp_px.clear()
            mcp_py.clear()
            mcp_pz.clear()
            mcp_id.clear()
            mcp_mother.clear()
            mcp_nuid.clear()
            mcp_pdg.clear()
            mcp_energy.clear()
            mcp_startx.clear()    
            mcp_starty.clear()
            mcp_startz.clear()
            mcp_endx.clear()
            mcp_endy.clear()
            mcp_endz.clear()
            nue.clear()
            nuID.clear()
            nuPDG.clear()
            nuvtxx.clear()
            nuvtxy.clear()
            nuvtxz.clear()
            nupx.clear()
            nupy.clear()
            nupz.clear()
            mode.clear()
            ccnc.clear()
            # setup vector to fill vector of vectors
            packetFrac=ROOT.std.vector("float")() 
            trackID=ROOT.std.vector("int")() 
            trackIndex=ROOT.std.vector("int")() 
            particleID=ROOT.std.vector("int")() 
            particleIndex=ROOT.std.vector("int")()
            interactionIndex=ROOT.std.vector("int")()
            pdgHit=ROOT.std.vector("int")()
            print('ev index of loop',ev_index ,len(events),end='\r')
            # Get event info for data
            event = f["charge/events/data"][ev_index]
            event_calib_prompt_hits=flow_out["charge/events/","charge/calib_prompt_hits", events["id"][ev_index]]
            event_calib_final_hits=flow_out["charge/events/","charge/calib_final_hits", events["id"][ev_index]]
            
            if (useData==False): 
                # find spillID to use for truth info
                spillArray=flow_out["charge/calib_prompt_hits","charge/packets","mc_truth/segments",event_calib_prompt_hits[0]["id"]]["event_id"][0][0][0]
                # find all truth info and fill it using a complicated vector 
                allTrajectories,allVertices=find_all_truth_in_spill(spillArray, flow_out)
                [nuID.push_back(int(i)) for i in allVertices[0]]
                [nue.push_back(i) for i in allVertices[1]]
                [nuPDG.push_back(int(i)) for i in allVertices[2]]
                [nuvtxx.push_back(i+trueXOffset) for i in allVertices[3]]
                [nuvtxy.push_back(i+trueYOffset) for i in allVertices[4]]
                [nuvtxz.push_back(i+trueZOffset) for i in allVertices[5]]
                [nupx.push_back(i) for i in allVertices[6]]
                [nupy.push_back(i) for i in allVertices[7]]
                [nupz.push_back(i) for i in allVertices[8]]
                [mode.push_back(i) for i in allVertices[9]]
                [ccnc.push_back(i) for i in allVertices[10]]
                [mcp_mother.push_back(int(i)) for i in allTrajectories[-1]]
                [mcp_startx.push_back(i+trueXOffset) for i in allTrajectories[0]]
                [mcp_starty.push_back(i+trueYOffset) for i in allTrajectories[1]]
                [mcp_startz.push_back(i+trueZOffset) for i in allTrajectories[2]]
                [mcp_endx.push_back(i+trueXOffset) for i in allTrajectories[3]]
                [mcp_endy.push_back(i+trueYOffset) for i in allTrajectories[4]]
                [mcp_endz.push_back(i+trueZOffset) for i in allTrajectories[5]]
                [mcp_px.push_back(i) for i in allTrajectories[6]]
                [mcp_py.push_back(i) for i in allTrajectories[7]]
                [mcp_pz.push_back(i) for i in allTrajectories[8]]
                [mcp_nuid.push_back(int(i)) for i in allTrajectories[-2]]
                [mcp_pdg.push_back(int(i)) for i in allTrajectories[-3]]
                [mcp_id.push_back(int(i)) for i in allTrajectories[-4]]
                [mcp_energy.push_back(i) for i in allTrajectories[-5]] 
            # fill event info
            eventID[0]           = event['id'] 
            event_start_t[0] =int(event["ts_start"])
            event_end_t[0]=int(event["ts_end"])
            event_unix_ts[0]=int(event["unix_ts"])
            # grab event hit list and variables
            hits_z=np.ma.getdata(event_calib_final_hits["z"][0])
            hits_y=np.ma.getdata(event_calib_final_hits["y"][0])
            hits_x=np.ma.getdata(event_calib_final_hits["x"][0])
            hits_Q=np.ma.getdata(event_calib_final_hits["Q"][0])
            hits_E=np.ma.getdata(event_calib_final_hits["E"][0])
            hits_ts=np.ma.getdata(event_calib_final_hits["ts_pps"][0])
            hits_id=np.ma.getdata(event_calib_final_hits["id"][0]) 
            hits_id_raw=np.ma.getdata(event_calib_final_hits["id"][0]) 
            # run over list of hits in the event
            hitID=0
            while hitID<len(hits_id): 
                packetFrac.clear()   
                trackID.clear()
                trackIndex.clear() 
                particleID.clear() 
                particleIndex.clear()
                interactionIndex.clear()
                interactionIndex.clear()
                pdgHit.clear()       
    

                hit_num=hits_id[hitID]
                contr_info=[]
                if (useData==False):              
                    # save 
                    contr_info=find_tracks_in_calib_hits(hit_num,flow_out)
                    [packetFrac.push_back(i) for i in contr_info[0]]
                    #[trackIndex.push_back(int(i)) for i in contr_info[1]]
                    [trackID.push_back(int(i)) for i in contr_info[2]]
                    [particleID.push_back(int(i)) for i in contr_info[4]]
                    #[particleIndex.push_back(int(i)) for i in contr_info[3]]
                    [pdgHit.push_back(int(i)) for i in contr_info[5]]
                # save hit information
                z.push_back(hits_z[hitID]+trueZOffset)
                y.push_back(hits_y[hitID]+trueYOffset)
                x.push_back(hits_x[hitID]+trueXOffset)
                charge.push_back(hits_Q[hitID])
                E.push_back(hits_E[hitID])
                ts.push_back(hits_ts[hitID])
                hit_packetFrac.push_back(packetFrac)
                hit_particleID.push_back(particleID)
                hit_particleIndex.push_back(particleIndex)
                hit_pdg.push_back(pdgHit)
                hit_interactionIndex.push_back(interactionIndex)
                hit_segmentIndex.push_back(trackIndex)
                hit_segmentID.push_back(trackID)

                hitID=hitID+1
            # grab event hit list and variables for prompt hits
            hits_z=np.ma.getdata(event_calib_prompt_hits["z"][0])
            hits_y=np.ma.getdata(event_calib_prompt_hits["y"][0])
            hits_x=np.ma.getdata(event_calib_prompt_hits["x"][0])
            hits_Q=np.ma.getdata(event_calib_prompt_hits["Q"][0])
            hits_E=np.ma.getdata(event_calib_prompt_hits["E"][0])
            hits_ts=np.ma.getdata(event_calib_prompt_hits["ts_pps"][0])
            hits_id=np.ma.getdata(event_calib_prompt_hits["id"][0]) 
            hits_id_raw=np.ma.getdata(event_calib_prompt_hits["id"][0]) 
            # run over list of hits in the event
            hitID=0
            while hitID<len(hits_id): 
                packetFrac.clear()   
                trackID.clear()
                trackIndex.clear() 
                particleID.clear() 
                particleIndex.clear()
                interactionIndex.clear()
                interactionIndex.clear()
                pdgHit.clear()       
    

                hit_num=hits_id[hitID]
                contr_info=[]
                if (useData==False):              
                    # save 
                    contr_info=find_tracks_in_calib_hits(hit_num,flow_out, typ="prompt")
                    [packetFrac.push_back(i) for i in contr_info[0]]
                    #[trackIndex.push_back(int(i)) for i in contr_info[1]]
                    [trackID.push_back(int(i)) for i in contr_info[2]]
                    [particleID.push_back(int(i)) for i in contr_info[4]]
                    #[particleIndex.push_back(int(i)) for i in contr_info[3]]
                    [pdgHit.push_back(int(i)) for i in contr_info[5]]
                # save hit information
                z_uncalib.push_back(hits_z[hitID]+trueZOffset)
                y_uncalib.push_back(hits_y[hitID]+trueYOffset)
                x_uncalib.push_back(hits_x[hitID]+trueXOffset)
                charge_uncalib.push_back(hits_Q[hitID])
                E_uncalib.push_back(hits_E[hitID])
                ts_uncalib.push_back(hits_ts[hitID])
                hit_packetFrac_uncalib.push_back(packetFrac)
                hit_particleID_uncalib.push_back(particleID)
                hit_particleIndex_uncalib.push_back(particleIndex)
                hit_pdg_uncalib.push_back(pdgHit)
                hit_interactionIndex_uncalib.push_back(interactionIndex)
                hit_segmentIndex_uncalib.push_back(trackIndex)
                hit_segmentID_uncalib.push_back(trackID)

                hitID=hitID+1
               
            output_tree.Fill()

        #end event loop

        #create indexx


        #####
        # Write output_tree
        # --------------------------------------------------------
        output_file.cd()
        output_tree.Write()


        #print('\n Data has been written to %s ' %(outputpath + '/' + outFileName))
    #end file loop
    #print("end of code")
    print("end of code")
def find_tracks_in_calib_hits(hit_num, flow_out, typ="final"):
    final_hit_backtrack=flow_out["charge/calib_"+str(typ)+"_hits","mc_truth/calib_"+str(typ)+"_hit_backtrack",hit_num][0]
    final_hit_backtrackFull=flow_out["charge/calib_"+str(typ)+"_hits","mc_truth/calib_"+str(typ)+"_hit_backtrack",hit_num]
    segIDsFromHits=final_hit_backtrack["segment_ids"][0]
    fracFromHits=final_hit_backtrack["fraction"][0]
    track_contr = []
    trackIndex_tot=[]
    segment_tot=[]
    pdg_tot=[] 
    particleIndex_tot=[]
    particle_tot=[]
    vertex_tot=[]
    vertexID_tot=[]
    pdg_id=[]
    total = 0.
    # Get fraction information and track information from hit
    i=0 
    while i<len(fracFromHits):
        fracs=fracFromHits[i]
        ids=segIDsFromHits[i]
        if (fracs==0 and ids==0):
            i=i+1 
            continue;
        total += fracs
        # get fraction information
        traj_index =0
        interaction_index=0
        contrib = fracs
        track_contr.append(contrib)

        # get track info
        seg=flow_out["mc_truth/segments/data"][ids]
        segID = seg["segment_id"]
        pdg = seg["pdg_id"]
        particleID = seg["traj_id"]
        vertexID=seg["vertex_id"]
        pdg_id.append(pdg)
        particleIndex_tot.append(traj_index)
        trackIndex_tot.append(-999)
        particle_tot.append(particleID)
        vertexID_tot.append(vertexID)
        vertex_tot.append(interaction_index)
        i=i+1
    if len(fracFromHits)<1:
        pdg_id.append(-999)
        particleIndex_tot.append(-999)
        trackIndex_tot.append(-999)
        particle_tot.append(-999)
        vertexID_tot.append(-999)
        vertex_tot.append(-999)
    if len(segIDsFromHits)<1:
        track_contr.append(-999)
    
     
    return [track_contr,trackIndex_tot,segment_tot,particleIndex_tot,particle_tot,pdg_id,vertexID_tot,vertex_tot]

def find_tracks_in_packet(hit_num, flow_out):
    # variables we wil need for later
    track_contr = []
    trackIndex_tot=[]
    segment_tot=[]
    pdg_tot=[] 
    particleIndex_tot=[]
    particle_tot=[]
    vertex_tot=[]
    vertexID_tot=[]
    pdg_id=[]
    total = 0.
    # Get fraction information and track information from hit
    trajFromHits=flow_out["charge/calib_prompt_hits","charge/packets","mc_truth/segments",hit_num][0][0]
    fracFromHits=flow_out["charge/calib_prompt_hits","charge/packets","mc_truth/packet_fraction",hit_num][0][0]
    for fracs in fracFromHits:
        total += fracs
        # get fraction information
        traj_index =0
        interaction_index=0
        contrib = fracs
        track_contr.append(contrib)
    for trajs in trajFromHits:
        # get track info
        seg = trajs["segment_id"]
        pdg = trajs["pdg_id"]
        particleID = trajs["traj_id"]
        vertexID=trajs["vertex_id"]
        pdg_id.append(pdg)
        particleIndex_tot.append(traj_index)
        trackIndex_tot.append(-999)
        particle_tot.append(particleID)
        vertexID_tot.append(vertexID)
        vertex_tot.append(interaction_index)
    if len(fracFromHits)<1:
        pdg_id.append(-999)
        particleIndex_tot.append(-999)
        trackIndex_tot.append(-999)
        particle_tot.append(-999)
        vertexID_tot.append(-999)
        vertex_tot.append(-999)
    if len(trajFromHits)<1:
        track_contr.append(-999)
    
     
    return [track_contr,trackIndex_tot,segment_tot,particleIndex_tot,particle_tot,pdg_id,vertexID_tot,vertex_tot]

def find_all_truth_in_spill(spillID, flow_out):
    trajStartX=[]
    trajStartY=[]
    trajStartZ=[]
    trajEndX=[]
    trajEndY=[]
    trajEndZ=[]
    trajPx=[]
    trajPy=[]
    trajPz=[]
    trajE=[]
    trajID=[]
    trajPDG=[]
    trajVertexID=[]
    trajParentID=[]
    traj_indicesArray = np.where(flow_out['mc_truth/trajectories/data']["event_id"] == spillID)[0] 
    # get all the mcparticle information
    for traj_indices in traj_indicesArray:
        trajStartX.append(flow_out["mc_truth/trajectories/data"][traj_indices]["xyz_start"][0])
        trajStartY.append(flow_out["mc_truth/trajectories/data"][traj_indices]["xyz_start"][1])
        trajStartZ.append(flow_out["mc_truth/trajectories/data"][traj_indices]["xyz_start"][2])
        trajEndX.append(flow_out["mc_truth/trajectories/data"][traj_indices]["xyz_end"][0])
        trajEndY.append(flow_out["mc_truth/trajectories/data"][traj_indices]["xyz_end"][1])
        trajEndZ.append(flow_out["mc_truth/trajectories/data"][traj_indices]["xyz_end"][2])
        trajID.append(flow_out["mc_truth/trajectories/data"][traj_indices]["traj_id"])
        trajPDG.append(flow_out["mc_truth/trajectories/data"][traj_indices]["pdg_id"])
        pdg=flow_out["mc_truth/trajectories/data"][traj_indices]["pdg_id"]
        px=flow_out["mc_truth/trajectories/data"][traj_indices]["pxyz_start"][0]*0.001
        py=flow_out["mc_truth/trajectories/data"][traj_indices]["pxyz_start"][1]*0.001
        pz=flow_out["mc_truth/trajectories/data"][traj_indices]["pxyz_start"][2]*0.001
        trajE.append(flow_out["mc_truth/trajectories/data"][traj_indices]["E_start"]*0.001)

        trajPx.append(px)
        trajPy.append(py)
        trajPz.append(pz)
        trajVertexID.append(flow_out["mc_truth/trajectories/data"][traj_indices]["vertex_id"])
        trajParentID.append(flow_out["mc_truth/trajectories/data"][traj_indices]["parent_id"])
    nuVertexIndex=[]
    for i in trajVertexID:
        if i<1000000:
            nuVertexIndex.append(int(i))
        else:
            a=str(i)
            a=int(a[0]+a[-5:])
            nuVertexIndex.append(a) 
    trajectories=[trajStartX,trajStartY,trajStartZ,trajEndX,trajEndY,trajEndZ,trajPx,trajPy,trajPz,trajE,trajID,trajPDG,nuVertexIndex,trajParentID]
    vertex_indicesArray = np.where(flow_out["/mc_truth/interactions/data"]["event_id"] == spillID)[0]
    # get all the neutrino information
    nuVertexID=[]
    nuVertexX=[]
    nuVertexY=[]
    nuVertexZ=[]
    nuVertexE=[]
    nuPDG=[]
    nuPx=[]
    nuPy=[]
    nuPz=[]
    nuCode=[]
    nuCC=[]
    nuVertexArray=[]
    for vertex_indices in vertex_indicesArray:
        nuVertexArray.append(flow_out["/mc_truth/interactions/data"][vertex_indices]["vertex_id"])
        nuVertexX.append(flow_out["/mc_truth/interactions/data"][vertex_indices]["x_vert"])
        nuVertexY.append(flow_out["/mc_truth/interactions/data"][vertex_indices]["y_vert"])
        nuVertexZ.append(flow_out["/mc_truth/interactions/data"][vertex_indices]["z_vert"])
        nuVertexE.append(flow_out["/mc_truth/interactions/data"][vertex_indices]["Enu"]*0.001)
        nuPDG.append(flow_out["/mc_truth/interactions/data"][vertex_indices]["nu_pdg"])
        nuPx.append(flow_out["/mc_truth/interactions/data"][vertex_indices]["nu_4mom"][0]*0.001)
        nuPy.append(flow_out["/mc_truth/interactions/data"][vertex_indices]["nu_4mom"][1]*0.01)
        nuPz.append(flow_out["/mc_truth/interactions/data"][vertex_indices]["nu_4mom"][2]*0.01)
        code,cc=get_nuance_code(vertex_indices,flow_out)
        nuCode.append(code)
        nuCC.append(cc)



    for i in nuVertexArray:
        if i<1000000:
            nuVertexID.append(int(i))
        else:
            a=str(i)
            a=int(a[0]+a[-5:])
            nuVertexID.append(a)
    vertices=[nuVertexID,nuVertexE,nuPDG,nuVertexX,nuVertexY,nuVertexZ,nuPx,nuPy,nuPz,nuCode,nuCC]
    return trajectories, vertices
def get_nuance_code(vertex_num,flow_out):
    # convert the neutrino information to nuance code from PandoraInterface code
    qe=flow_out["mc_truth/interactions/data"][vertex_num]["isQES"]
    nccc=flow_out["mc_truth/interactions/data"][vertex_num]["isCC"]
    mec=flow_out["mc_truth/interactions/data"][vertex_num]["isMEC"]
    res=flow_out["mc_truth/interactions/data"][vertex_num]["isRES"]
    dis=flow_out["mc_truth/interactions/data"][vertex_num]["isDIS"]
    coh=flow_out["mc_truth/interactions/data"][vertex_num]["isCOH"]
    code=1000;
    cc=1
    if nccc==1:
        cc=0
    if (qe):
        code=0

    if (dis):
        code=2
    if (res):
        code=1
    if (coh):
        code=3
        if (qe):
            code=4;
    if (mec):
        code=10
    return int(code),int(cc)
if __name__ == "__main__":
    main()

