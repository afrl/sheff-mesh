import mesh

m = mesh.Mesh()
m.nodes().add([0, 0, 0])

for n in m.nodes():
    print('Node ', n.local_index())
    print('  Local index: ', n.local_index())
    print('  Global index: ', n.get_global_index())
    print('  Coordinate: ', list(n.get_point()))
    print('  Partition: ', n.get_parent_partition())

for n in m.nodes():
    n.get_point()[0] = 1.0

print('\nAfter modication of coordinates')
for n in m.nodes():
    print('Node ', n.local_index())
    print('  Local index: ', n.local_index())
    print('  Global index: ', n.get_global_index())
    print('  Coordinate: ', list(n.get_point()))
    print('  Partition: ', n.get_parent_partition())