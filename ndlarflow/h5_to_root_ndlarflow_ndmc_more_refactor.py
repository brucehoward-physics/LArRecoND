import numpy as np
import awkward as awk

import uproot as ur
import h5py as h5

from h5flow.data import dereference
import h5flow
import os
import sys

# Refactored version of h5 to ROOT conversion, Bruce Howard - 2025

# NOTE for more on uproot TTree writing see the below or the uproot documentation
# see https://stackoverflow.com/questions/72187937/writing-trees-number-of-baskets-and-compression-uproot

# Main function with command line settable params
def main(argv=None):
    fileNames=[]
    useData=False
    useFinalHits=False
    overrideOutname=1
    outname=''

    MeV2GeV=0.001
    trueXOffset=0 # Offsets if geometry changes
    trueYOffset=0#42+268
    trueZOffset=0#-1300

    if len(sys.argv==1):
        print('Must at least pass a file location/name to be converted')
    if len(sys.argv)>1:
        if sys.argv[1]!=None:
            fileList=str(sys.argv[1])
            fileNames=fileList.split(',')
        if len(sys.argv)>2 and sys.argv[2]!=None:
            if int(sys.argv[2])==1:
                useData=True
        if len(sys.argv)>3 and sys.argv[3]!=None:
            if int(sys.argv[3])==1:
                useFinalHits=True
        if len(sys.argv)>4 and sys.argv[4]!=None:
            outname=str(sys.argv[4])
            overrideOutname=0

    MaxArrayDepth=int(10000)
    isWritten=False

    promptKey='prompt'
    if useFinalHits==True:
        promptKey='final'

    if overrideOutname==1:
        outname = fileNames[0].split('/')[-1]+'_hits.root'

    for fileIdx in range(len(fileNames)):
        print('Processing file',fileIdx,'of',len(fileNames))
        fileName=fileNames[fileIdx]

        f = h5.File(fileName)
        events=f['charge/events/data']
        flow_out=h5flow.data.H5FlowDataManager(fileName,"r")

        eventsToRun=len(events)

        for ievt in range(eventsToRun):
            if ievt%10==0:
                print('Currently on',ievt,'of',eventsToRun)
            event = events[ievt]
            event_calib_prompt_hits=flow_out["charge/events/","charge/calib_"+promptKey+"_hits", events["id"][ievt]]

            if len(event_calib_prompt_hits)==0:
                print('This event seems empty in the hits array, skipping')
                continue

            # Start with the non-spill info, this is all ~like the current form
            # note, not sure why this is repeated in the main python code...
            runID = np.array( [0], dtype='int32' )
            subrunID = np.array( [0], dtype='int32' )
            eventID = np.array( [event['id']], dtype='int32' )
            event_start_t = np.array( [event['ts_start']], dtype='int32' )
            event_end_t = np.array( [event['ts_end']], dtype='int32' )
            event_unix_ts = np.array( [event['unix_ts']], dtype='int32' )

            # Removing duplicate hits_id instantiation and getting rid of hits_id_raw which is unused
            #######################################
            hits_z = (np.ma.getdata(event_calib_prompt_hits["z"][0])+trueZOffset).astype('float32')
            hits_y = ( np.ma.getdata(event_calib_prompt_hits["y"][0])+trueYOffset ).astype('float32')
            hits_x = ( np.ma.getdata(event_calib_prompt_hits["x"][0])+trueXOffset ).astype('float32')
            hits_Q = ( np.ma.getdata(event_calib_prompt_hits["Q"][0]) ).astype('float32')
            hits_E = ( np.ma.getdata(event_calib_prompt_hits["E"][0]) ).astype('float32')
            hits_ts = ( np.ma.getdata(event_calib_prompt_hits["ts_pps"][0]) ).astype('float32')
            hits_ids = np.ma.getdata(event_calib_prompt_hits["id"][0])

            # "uncalib" -- same for now, not even used in LArPandora... why save? Make it optional to use the prompt or final hits and leave at that?
            #######################################

            if len(hits_ids)==0:
                print('This event has no hit IDs, skipping')
                continue

            if useData==False:
                # Truth-level info for hits
                #######################################
                trajFromHits=flow_out["charge/calib_prompt_hits","charge/packets","mc_truth/segments",hits_ids[:]][:,0]
                fracFromHits=flow_out["charge/calib_prompt_hits","charge/packets","mc_truth/packet_fraction",hits_ids[:]][:,0]

                matches = trajFromHits['segment_id'].count(axis=1).astype('uint16')

                packetFrac = np.array( awk.flatten( awk.flatten( awk.Array( fracFromHits['fraction'].astype('float32') ) ) ) )
                packetFrac = packetFrac[np.where(packetFrac!=0)]

                trajFromHits=trajFromHits.data[~trajFromHits['segment_id'].mask]
                pdgHit = trajFromHits['pdg_id'].astype('int32')
                trackID = trajFromHits['segment_id'].astype('int32')
                particleID = trajFromHits['file_traj_id'].astype('int64')
                particleIDLocal = trajFromHits['traj_id'].astype('int64')
                interactionIndex = trajFromHits['vertex_id'].astype('int64')

                # Truth-level info for the spill
                #######################################
                spillID=flow_out["charge/calib_prompt_hits","charge/packets","mc_truth/segments",hits_ids]["event_id"][0][0][0]
                # Trajectories
                traj_indicesArray = np.where(flow_out['mc_truth/trajectories/data']["event_id"] == spillID)[0]
                traj = flow_out["mc_truth/trajectories/data"][traj_indicesArray]
                trajStartX = (traj['xyz_start'][:,0]).astype('float32')
                trajStartY = (traj['xyz_start'][:,1]).astype('float32')
                trajStartZ = (traj['xyz_start'][:,2]).astype('float32')
                trajEndX = (traj['xyz_end'][:,0]).astype('float32')
                trajEndY = (traj['xyz_end'][:,1]).astype('float32')
                trajEndZ = (traj['xyz_end'][:,2]).astype('float32')
                trajLength = (traj['dist_travel']).astype('float32')
                trajTStart = (traj['t_start']).astype('double')
                trajTEnd = (traj['t_end']).astype('double')
                trajID = (traj['file_traj_id']).astype('int64')
                trajIDLocal = (traj['traj_id']).astype('int64')
                trajPDG = (traj['pdg_id']).astype('int32')
                trajE = (traj['E_start']*MeV2GeV).astype('float32')
                trajPx = (traj['pxyz_start'][:,0]*MeV2GeV).astype('float32')
                trajPy = (traj['pxyz_start'][:,1]*MeV2GeV).astype('float32')
                trajPz = (traj['pxyz_start'][:,2]*MeV2GeV).astype('float32')
                trajVertexID = (traj['vertex_id']).astype('int64')
                trajParentID = (traj['parent_id']).astype('int64')
                # Vertices
                vertex_indicesArray = np.where(flow_out["/mc_truth/interactions/data"]["event_id"] == spillID)[0]
                vtx = flow_out["/mc_truth/interactions/data"][vertex_indicesArray]
                nu_vtx_id = (vtx['vertex_id']).astype('int64')
                nu_vtx_x = (vtx['x_vert']).astype('float32')
                nu_vtx_y = (vtx['y_vert']).astype('float32')
                nu_vtx_z = (vtx['z_vert']).astype('float32')
                nu_vtx_E = (vtx['Enu']*MeV2GeV).astype('float32')
                nu_spill_t = (vtx['t_event']).astype('float64')
                nu_pdg = (vtx['nu_pdg']).astype('int32')
                nu_px = (vtx['nu_4mom'][:,0]*MeV2GeV).astype('float32')
                nu_py = (vtx['nu_4mom'][:,1]*MeV2GeV).astype('float32')
                nu_pz = (vtx['nu_4mom'][:,2]*MeV2GeV).astype('float32')
                # Little bit of gymnastics here
                ccnc = vtx['isCC']
                nu_iscc = np.invert(ccnc).astype('int32')
                # And more gymnastics here
                codes = 1000*np.ones(len(nu_vtx_id),dtype='int32')
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
                nu_code = codes

            ## Rebuild now with all the individual types
            event_dict = { 'run':runID, 'subrun':subrunID, 'event':eventID, 'unix_ts':event_unix_ts, 'event_start_t':event_start_t, 'event_end_t':event_end_t }

            if useData==False:
                other_dict = {  'x':hits_x, 'y':hits_y, 'z':hits_z, 'ts':hits_ts, 'charge':hits_Q, 'E':hits_E, 'matches':matches,\
                                'mcp_energy':trajE, 'mcp_pdg':trajPDG, 'mcp_nuid':trajVertexID, 'mcp_vertex_id':trajVertexID,\
                                'mcp_idLocal':trajIDLocal, 'mcp_id':trajID, 'mcp_px':trajPx, 'mcp_py':trajPy, 'mcp_pz':trajPz,\
                                'mcp_mother':trajParentID, 'mcp_startx':trajStartX, 'mcp_starty':trajStartY, 'mcp_startz':trajStartZ,\
                                'mcp_endx':trajEndX, 'mcp_endy':trajEndY, 'mcp_endz':trajEndZ, 'mcp_length':trajLength, 'mcp_tstart':trajTStart, 'mcp_tend':trajTEnd,\
                                'nuID':nu_vtx_id, 'vertex_id':nu_vtx_id, 'nue':nu_vtx_E, 'nuspillt':nu_spill_t, 'nuPDG':nu_pdg,\
                                'nupx':nu_px, 'nupy':nu_py, 'nupz':nu_pz, 'nuvtxx':nu_vtx_x, 'nuvtxy':nu_vtx_y,\
                                'nuvtxz':nu_vtx_z, 'mode':nu_code, 'ccnc':nu_iscc,\
                                'hit_packetFrac':packetFrac, 'hit_particleID':particleID, 'hit_particleIDLocal':particleIDLocal,\
                                'hit_pdg':pdgHit, 'hit_vertexID':interactionIndex, 'hit_segmentID':trackID }
            else:
                other_dict = {  'x':hits_x, 'y':hits_y, 'z':hits_z, 'ts':hits_ts, 'charge':hits_Q, 'E':hits_E }

            max_entries=0
            for key in other_dict.keys():
                if len(other_dict[key]) > max_entries:
                    max_entries = len(other_dict[key])

            nSubEvents = int(max_entries/MaxArrayDepth)+1
            for idxSubEvent in range(nSubEvents):
                first = MaxArrayDepth*idxSubEvent
                last = MaxArrayDepth*(idxSubEvent+1)
                event_dict['subevent'] = np.array([idxSubEvent], dtype='int32')
                for key in other_dict.keys():
                    event_dict[key] = awk.values_astype(awk.Array([other_dict[key][first:last]]),other_dict[key].dtype)
                if isWritten==False:
                    fout = ur.recreate(outname)
                    fout['events'] = event_dict
                    isWritten=True
                else:
                    fout['events'].extend(event_dict)

            del packetFrac
            del particleID
            del particleIDLocal
            del pdgHit
            del interactionIndex
            del trackID

        fout.close()
        print('end of code')

if __name__=="__main__":
    main()
