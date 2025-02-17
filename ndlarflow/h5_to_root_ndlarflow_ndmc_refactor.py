import numpy as np
import awkward as awk

import uproot as ur
import h5py as h5

from h5flow.data import dereference
import h5flow
import os
import sys

# NOTE
# see https://stackoverflow.com/questions/72187937/writing-trees-number-of-baskets-and-compression-uproot

## ONLY CHANGE THESE. MAKE THESE CMD LINE PARAMS

fileName='/exp/dune/data/users/howard/MicroProdN3p1_NDLAr_2E18_FHC.flow.nu.0000001.FLOW.hdf5'
eventsToRun=-1

### ------------------ DO NOT CHANGE BELOW HERE

f = h5.File(fileName)
outname = fileName.split('/')[-1]+'_hits.root'
events=f['charge/events/data']
flow_out=h5flow.data.H5FlowDataManager(fileName,"r")

if eventsToRun < 0:
    eventsToRun = len(events)

MeV2GeV=0.001
trueXOffset=0 # Offsets if geometry changes
trueYOffset=0#42+268
trueZOffset=0#-1300

maxhits=0

for ievt in range(eventsToRun):
    if ievt%10==0:
        print('Currently on',ievt,'of',eventsToRun)
    event = events[ievt]
    event_calib_prompt_hits=flow_out["charge/events/","charge/calib_prompt_hits", events["id"][ievt]]

    # Start with the non-spill info, this is all ~like the current form
    # note, not sure why this is repeated in the main python code...
    runID = np.array( [0], dtype='int' )
    subrunID = np.array( [0], dtype='int' )
    eventID = np.array( [event['id']], dtype='int' )
    event_start_t = np.array( [event['ts_start']], dtype=np.int32 )
    event_end_t = np.array( [event['ts_end']], dtype=np.int32 )
    event_unix_ts = np.array( [event['unix_ts']], dtype=np.int32 )

    # Removing duplicate hits_id instantiation and getting rid of hits_id_raw which is unused
    #######################################
    hits_z=awk.Array([ awk.Array( np.ma.getdata(event_calib_prompt_hits["z"][0])+trueZOffset ) ])
    hits_y=awk.Array([ awk.Array( np.ma.getdata(event_calib_prompt_hits["y"][0])+trueYOffset ) ])
    hits_x=awk.Array([ awk.Array( np.ma.getdata(event_calib_prompt_hits["x"][0])+trueXOffset ) ])
    hits_Q=awk.Array([ awk.Array( np.ma.getdata(event_calib_prompt_hits["Q"][0]) ) ])
    hits_E=awk.Array([ awk.Array( np.ma.getdata(event_calib_prompt_hits["E"][0]) ) ])
    hits_ts=awk.Array([ awk.Array( np.ma.getdata(event_calib_prompt_hits["ts_pps"][0]) ) ])
    hits_id=awk.Array([ awk.Array( np.ma.getdata(event_calib_prompt_hits["id"][0]) ) ])

    if len(hits_id[0]) > maxhits:
        maxhits = len(hits_id[0])

    # "uncalib" -- same for now
    #######################################
    hits_z_uncalib=awk.Array([ awk.Array( np.ma.getdata(event_calib_prompt_hits["z"][0])+trueZOffset ) ])
    hits_y_uncalib=awk.Array([ awk.Array( np.ma.getdata(event_calib_prompt_hits["y"][0])+trueYOffset ) ])
    hits_x_uncalib=awk.Array([ awk.Array( np.ma.getdata(event_calib_prompt_hits["x"][0])+trueXOffset ) ])
    hits_Q_uncalib=awk.Array([ awk.Array( np.ma.getdata(event_calib_prompt_hits["Q"][0]) ) ])
    hits_E_uncalib=awk.Array([ awk.Array( np.ma.getdata(event_calib_prompt_hits["E"][0]) ) ])
    hits_ts_uncalib=awk.Array([ awk.Array( np.ma.getdata(event_calib_prompt_hits["ts_pps"][0]) ) ])
    hits_id_uncalib=awk.Array([ awk.Array( np.ma.getdata(event_calib_prompt_hits["id"][0]) ) ])

    # Truth-level info for hits
    #######################################
    trajFromHits=flow_out["charge/calib_prompt_hits","charge/packets","mc_truth/segments",hits_id[0][:]][:,0]
    fracFromHits=flow_out["charge/calib_prompt_hits","charge/packets","mc_truth/packet_fraction",hits_id[0][:]][:,0]

    matches = awk.Array([ trajFromHits['segment_id'].count(axis=1) ])

    packetFrac = np.array( awk.flatten( awk.flatten( awk.Array( fracFromHits['fraction'].astype('float') ) ) ) )
    packetFrac = awk.Array([ packetFrac[ awk.where(packetFrac!=0) ] ])

    trajFromHits=trajFromHits.data[~trajFromHits['segment_id'].mask]
    pdgHit = awk.Array([ trajFromHits['pdg_id'].astype('int') ])
    trackID = awk.Array([ trajFromHits['segment_id'].astype('int') ])
    particleID = awk.Array([ trajFromHits['file_traj_id'].astype('int') ])
    particleIDLocal = awk.Array([ trajFromHits['traj_id'].astype('int') ])
    interactionIndex = awk.Array([ trajFromHits['vertex_id'].astype('int') ])

    # Truth-level info for the spill
    #######################################
    spillID=flow_out["charge/calib_prompt_hits","charge/packets","mc_truth/segments",hits_id[0]]["event_id"][0][0][0]
    # Trajectories
    traj_indicesArray = np.where(flow_out['mc_truth/trajectories/data']["event_id"] == spillID)[0]
    traj = flow_out["mc_truth/trajectories/data"][traj_indicesArray]
    trajStartX = awk.Array([ awk.Array(traj['xyz_start'][:,0]) ])
    trajStartY = awk.Array([ awk.Array(traj['xyz_start'][:,1]) ])
    trajStartZ = awk.Array([ awk.Array(traj['xyz_start'][:,2]) ])
    trajEndX = awk.Array([ awk.Array(traj['xyz_end'][:,0]) ])
    trajEndY = awk.Array([ awk.Array(traj['xyz_end'][:,1]) ])
    trajEndZ = awk.Array([ awk.Array(traj['xyz_end'][:,2]) ])
    trajLength = awk.Array([ awk.Array(traj['dist_travel']) ])
    trajTStart = awk.Array([ awk.Array(traj['t_start']) ])
    trajTEnd = awk.Array([ awk.Array(traj['t_end']) ])
    trajID = awk.Array([ awk.Array(traj['file_traj_id']) ])
    trajIDLocal = awk.Array([ awk.Array(traj['traj_id']) ])
    trajPDG = awk.Array([ awk.Array(traj['pdg_id']) ])
    trajE = awk.Array([ awk.Array(traj['E_start']*MeV2GeV) ])
    trajPx = awk.Array([ awk.Array(traj['pxyz_start'][:,0]*MeV2GeV) ])
    trajPy = awk.Array([ awk.Array(traj['pxyz_start'][:,1]*MeV2GeV) ])
    trajPz = awk.Array([ awk.Array(traj['pxyz_start'][:,2]*MeV2GeV) ])
    trajVertexID = awk.Array([ awk.Array(traj['vertex_id']) ])
    trajParentID = awk.Array([ awk.Array(traj['parent_id']) ])
    # Vertices
    vertex_indicesArray = np.where(flow_out["/mc_truth/interactions/data"]["event_id"] == spillID)[0]
    vtx = flow_out["/mc_truth/interactions/data"][vertex_indicesArray]
    nu_vtx_id = awk.Array([ awk.Array(vtx['vertex_id']) ])
    nu_vtx_x = awk.Array([ awk.Array(vtx['x_vert']) ])
    nu_vtx_y = awk.Array([ awk.Array(vtx['y_vert']) ])
    nu_vtx_z = awk.Array([ awk.Array(vtx['z_vert']) ])
    nu_vtx_E = awk.Array([ awk.Array(vtx['Enu']*MeV2GeV) ])
    nu_spill_t = awk.Array([ awk.Array(vtx['t_event']) ])
    nu_pdg = awk.Array([ awk.Array(vtx['nu_pdg']) ])
    nu_px = awk.Array([ awk.Array(vtx['nu_4mom'][:,0]*MeV2GeV) ])
    nu_py = awk.Array([ awk.Array(vtx['nu_4mom'][:,1]*MeV2GeV) ])
    nu_pz = awk.Array([ awk.Array(vtx['nu_4mom'][:,2]*MeV2GeV) ])
    # Little bit of gymnastics here
    ccnc = vtx['isCC']
    nu_iscc = awk.Array([ awk.Array(np.invert(ccnc).astype('int')) ])
    # And more gymnastics here
    codes = 1000*np.ones(len(nu_vtx_id[0]),dtype='int')
    idxQE = np.where(vtx['isQES']==True)
    idxRES = np.where(vtx['isRES']==True)
    idxDIS = np.where(vtx['isDIS']==True)
    idxMEC = np.where(vtx['isMEC']==True)
    idxCOH = np.where(vtx['isCOH']==True)
    idxCOHQE = np.where((vtx['isCOH']==True) & (vtx['isQES']==True))
    codes[idxQE] = 0
    codes[idxRES] = 1
    codes[idxDIS] = 2
    codes[idxCOH] = 3
    codes[idxCOHQE] = 4
    codes[idxMEC] = 10
    nu_code = awk.Array([ awk.Array(codes) ])

    # Putting this into a dict to write to TTree:
    print(len(matches),len(hits_x))
    hits_dict = awk.zip({ 'x' : hits_x, 'y' : hits_y, 'z' : hits_z, 'ts' : hits_ts, 'charge' : hits_Q, 'E' : hits_E,\
                          'x_uncalib' : hits_x, 'y_uncalib' : hits_y, 'z_uncalib' : hits_z, 'ts_uncalib' : hits_ts,\
                          'charge_uncalib' : hits_Q, 'E_uncalib' : hits_E, 'matches' : matches })
    mcp_dict = awk.zip({ 'mcp_energy' : trajE, 'mcp_pdg' : trajPDG, 'mcp_nuid' : trajVertexID, 'mcp_vertex_id' : trajVertexID,\
                         'mcp_idLocal' : trajIDLocal, 'mcp_id' : trajID, 'mcp_px' : trajPx, 'mcp_py' : trajPy, 'mcp_pz' : trajPz,\
                         'mcp_mother' : trajParentID, 'mcp_startx' : trajStartX, 'mcp_starty' : trajStartY, 'mcp_startz' : trajStartZ,\
                         'mcp_endx' : trajEndX, 'mcp_endy' : trajEndY, 'mcp_endz' : trajEndZ, 'mcp_length' : trajLength,\
                         'mcp_tstart' : trajTStart, 'mcp_tend' : trajTEnd })
    nu_dict = awk.zip({ 'nuID' : nu_vtx_id, 'vertex_id' : nu_vtx_id, 'nue' : nu_vtx_E, 'nuspillt' : nu_spill_t, 'nuPDG' : nu_pdg,\
                        'nupx' : nu_px, 'nupy' : nu_py, 'nupz' : nu_pz, 'nuvtxx' : nu_vtx_x, 'nuvtxy' : nu_vtx_y,\
                        'nuvtxz' : nu_vtx_z, 'mode' : nu_code, 'ccnc' : nu_iscc })
    # Track index and particle index are left unfilled I think...
    event_dict = { 'run' : runID, 'subrun' : subrunID, 'event' : eventID, 'unix_ts' : event_unix_ts,\
                   'event_start_t' : event_start_t, 'event_end_t' : event_end_t,\
                   'hits' : hits_dict, 'mcparticles' : mcp_dict, 'neutrinos' : nu_dict,\
		   'hit_packetFrac' : packetFrac, 'hit_particleID' : particleID, 'hit_particleIDLocal' : particleIDLocal,\
		   'hit_pdg' : pdgHit, 'hit_vertexID' : interactionIndex, 'hit_segmentID' : trackID,\
		   'hit_packetFrac_uncalib' : packetFrac, 'hit_particleID_uncalib' : particleID, 'hit_particleIDLocal_uncalib' : particleIDLocal,\
		   'hit_pdg_uncalib' : pdgHit, 'hit_vertexID_uncalib' : interactionIndex, 'hit_segmentID_uncalib' : trackID }

    if ievt == 0:
        fout = ur.recreate(outname)
        fout['events'] = event_dict
    else:
        fout['events'].extend(event_dict)

    # Delete stuff that is hogging memory
    del packetFrac
    del particleID
    del particleIDLocal
    del pdgHit
    del interactionIndex
    del trackID

fout['maxhits'] = {'maxhits' : [maxhits]}
fout.close()
