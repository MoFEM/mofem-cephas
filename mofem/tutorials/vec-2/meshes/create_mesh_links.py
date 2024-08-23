cubit.cmd('reset')
cubit.cmd('set developer commands on')
cubit.cmd('set duplicate block elements on')

one_master_node = False

cubit.cmd('brick x 0.5 y 0.05 z 0.02')

cubit.cmd('block 1 add volume 1')
cubit.cmd('block 1 name "BODY_FORCE"')
cubit.cmd('block 1 attribute count 3')
cubit.cmd('block 1 attribute index 1 0')
cubit.cmd('block 1 attribute index 2 -1 ')
cubit.cmd('block 1 attribute index 2 0')

cubit.cmd('block 3 add surface 4 ')
cubit.cmd('block 3 name "FIX_ALL"')

cubit.cmd('block 4 add volume 1')
cubit.cmd('block 4 name "MAT_ELASTIC"')
cubit.cmd('block 4 attribute count 2')
cubit.cmd('block 4 attribute index 1 100')
cubit.cmd('block 4 attribute index 2 0.3')

cubit.cmd('webcut volume 1 with plane xplane offset 0.05')
cubit.cmd('webcut volume 2 with plane xplane offset 0')

cubit.cmd('delete volume 2')

cubit.cmd('volume all scheme tetmesh')
cubit.cmd('surface 22  size auto factor 9')
cubit.cmd('mesh surface 22 ')
cubit.cmd('volume all size auto factor 8')

cubit.cmd('#copy mesh surface 22  onto  surface 7  source curve 38  source vertex 22  target curve 15  target vertex 12   mirror ')
cubit.cmd('#mesh surface 7')

cubit.cmd('surface 7 scheme copy source surface 22 source vertex 22  target vertex 12 source curve 38  target curve 15   ')
cubit.cmd('mesh surface 7')
cubit.cmd('mesh volume all')


node_list1 = cubit.parse_cubit_list("node", "all in surface 22")	
node_list2 = cubit.parse_cubit_list("node", "all in surface 7")	

if len(node_list1) != len(node_list2):
    print("ERROR: node list size in surface 7 and surface 22 are different")

if one_master_node == True:
    node_list2 = node_list2[:1]
    # cubit.cmd('delete volume 1')


edges_ids = []
curve_ids = []

for id1 in node_list1:
  shortest_edge = 1e16
  shortest_id2 = None
  for id2 in node_list2:
    cubit.cmd('create edge node {} {}'.format(id1, id2))
    edge_id = cubit.get_last_id("edge")
    dist = cubit.get_mesh_edge_length(edge_id)
    cubit.cmd('delete edge {}'.format(edge_id))
    if dist < shortest_edge:
      shortest_id2 = id2
      shortest_edge = dist
  if shortest_id2 is not None:      
    cubit.cmd('create curve location at node {} location at node {} '.format(id1, shortest_id2))	
    curve_id = cubit.get_last_id("curve")
    curve_ids.append(curve_id)
    #cubit.cmd('create edge node {} {} owner curve {}'.format(id1, shortest_id2, curve_id))
    # cubit.cmd('create edge node {} {}'.format(id1, shortest_id2))
    # cubit.cmd('block 4 joint node {} spider node {} element type beam'.format(id1, shortest_id2))

    edge_id = cubit.get_last_id("edge")
    edges_ids.append(edge_id)


print(curve_ids)
print(edges_ids)

print("we are done")
cubit.cmd('merge all')

curve_list  = ' '.join(map(str, curve_ids))

cubit.cmd('curve {} interval 1'.format(curve_list))
cubit.cmd('mesh curve {}'.format(curve_list))
cubit.cmd('')

#cubit.cmd('disassociate mesh from volume all ')
#cubit.cmd('disassociate mesh from curve all ')
cubit.cmd('equivalence node all tolerance 0.001 force')
#cubit.cmd('delete volume all  ')
#cubit.cmd('delete curve all  ')

edge_list2 = cubit.parse_cubit_list("edge", "all in curve {} ".format(curve_list))
edge_list  = ' '.join(map(str, edges_ids))
edge_list  = ' '.join(map(str, edge_list2))

cubit.cmd('block 4 add edge {}'.format(edge_list))

cubit.cmd('block 4 name "MPC_COUPLING_LINKS4"')

#renumber Edge all start_id 99999
#renumber Edge all start_id 1


cubit.cmd('save as "/home/my_user/Desktop/beam3D_MPCs.cub" overwrite')

if one_master_node == True:
    cubit.cmd('delete volume 1')
    cubit.cmd('create force on vertex 26 vector 1 0 0 0 0 0')
    cubit.cmd('block 5 add node {}'.format(22))
    cubit.cmd('block 5 add vertex {}'.format(26))
    cubit.cmd('block 5 name "MPC_COUPLING_MASTER4"')
    cubit.cmd('save as "/home/my_user/Desktop/beam3D_MPCs_MASTER_SINGLE.cub" overwrite')


#cubit.cmd('block 5 add node all in surface {}'.format(22))
cubit.cmd('block 5 add surface {}'.format(22))
cubit.cmd('block 5 name "MPC_COUPLING_MASTER4"')

#cubit.cmd('block 6 add surface {}'.format(7))
#cubit.cmd('block 6 name "MPC_COUPLING_SLAVE4"')


cubit.cmd('save as "/home/my_user/Desktop/beam3D_MPCs_MASTER.cub" overwrite')
