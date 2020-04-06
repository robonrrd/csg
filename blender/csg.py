bl_info = {
    "name": "CSG",
    "category": "Object",
}


import bpy
import bmesh
import pycsg

"""
Utility Functions
"""
def is_mesh(obj):
    if obj is None:
        return False
    if obj.type in ['MESH']:
        return True
    return False


"""
Blender addon functionality
"""

class CSG(bpy.types.Operator):
    """Constructive solid geometry""" # blender will use this as a tooltip for menu items and buttons
    bl_idname = "object.csg" # unique identifier for buttons and menu items to reference.
    bl_label = "CSG"         # display name in the interface.
    bl_options = {'REGISTER', 'UNDO'}  # enable undo for the operator.

    def execute(self, context):
        scene = context.scene

        # CSG classes
        clay = None
        knife = None

        # Get selected blender objects.
        # There doesn't appear to be any way to get the blender
        # selection in order it was selected (!!!), so we need to
        # work around it.
        # We assume that the current active selection is the
        # knife (i.e. the clay is selected first)

        clay_obj = None
        knife_obj = None
        selected = bpy.context.selected_objects
        num_selected = len(selected)
        if num_selected != 2:
            raise Error("CSG need exactly two meshes selected (clay first, knife second)")

        if bpy.context.active_object == selected[1]:
            clay_obj = selected[0]
        else:
            clay_obj = selected[1]
        knife_obj = bpy.context.active_object

        if clay_obj.type in ['MESH']:
            print("Using",clay_obj.name,"as the clay")
        else:
            raise Error(clay_obj.name," is not a mesh")
        if knife_obj.type in ['MESH']:
            print("Using",knife_obj.name,"as the knife")
        else:
            raise Error(knife_obj.name," is not a mesh")

        clay = build_trimesh_from_obj(clay_obj)
        knife = build_trimesh_from_obj(knife_obj)

        if clay is not None and knife is not None:
            engine = pycsg.CSGEngine(clay, knife)
            A = pycsg.TriMesh()
            B = pycsg.TriMesh()
            engine.construct(pycsg.kDifference, True, A, B)
            build_blender_mesh_from_trimesh(A, "above")
            build_blender_mesh_from_trimesh(B, "below")

        return {'FINISHED'}            # this lets blender know the operator finished successfully.


def build_trimesh_from_obj(in_obj):
    if not is_mesh(in_obj):
        print("Passed in object is not a mesh. Failing.")
        return None

    faces = []
    for poly in in_obj.data.polygons:
        print("Polygon index: %d, length: %d" % (poly.index, poly.loop_total))
        for loop_index in range(poly.loop_start, poly.loop_start + poly.loop_total):
            faces.append(in_obj.data.loops[loop_index].vertex_index)

    vertices = []
    for vert in in_obj.data.vertices:
        pos = in_obj.matrix_world * vert.co
        vertices.extend(pos)

    #print("vertices", vertices)
    #print("faces", faces)
    mm = pycsg.TriMesh(vertices, faces)

    return mm


def build_blender_mesh_from_trimesh(trimesh, name):
    tm_verts = trimesh.mesh_vertices()
    num_verts = len(tm_verts)//3

    tm_faces = trimesh.mesh_faces()
    num_faces = len(tm_faces)//3

    print( "# vertices =", num_verts )
    print("# faces", num_faces)

    if num_verts == 0 or num_faces == 0:
        print("build_blender_mesh_from_trimesh: Empty mesh")
        return

    verts_data = []
    for i in range(num_verts):
        verts_data.append( (tm_verts[3*i], tm_verts[3*i+1], tm_verts[3*i+2]) )
    faces_data = []
    for i in range(num_faces):
        faces_data.append( (tm_faces[3*i], tm_faces[3*i+1], tm_faces[3*i+2]) )

    #print( "verts", verts_data )
    #print( "faces", faces_data )

    # create new mesh structure
    mesh = bpy.data.meshes.new(str(name+"_mesh"))
    mesh.from_pydata(verts_data, [], faces_data)
    mesh.update()

    new_object = bpy.data.objects.new(str(name+"_object"), mesh)
    new_object.data = mesh

    scene = bpy.context.scene
    scene.objects.link(new_object)
    scene.objects.active = new_object
    new_object.select = True



    '''
    # works
    vertsData = [(-1.509843, -1.525169, 0.433532),
                (-1.509843, 0.474831, 0.433532),
                (-1.509843, -1.525169, -1.566468),
                (-1.509843, 0.474831, -1.566468),
                (0.490157, -1.525169, 0.433532),
                (0.490157, 0.474831, 0.433532),
                (0.490157, -1.525169, -1.566468),
                (0.490157, 0.474831, -1.566468)]
    facesData = [(3, 2, 0), # 4 3 1
             (7, 6, 2), # 8 7 3
             (5, 3, 6), # 6 5 7
             (1, 0, 4), # 2 1 5
             (2, 6, 4), # 3 7 5
             (7, 3, 1), # 8 4 2
             (1, 3, 0), # 2 4 1
             (3, 7, 2), # 4 8 3
             (7, 5, 6), # 8 6 7
             (5, 1, 4), # 6 2 5
             (0, 2, 4), # 1 3 5
             (5, 7, 1)] # 6 8 2

    # create new mesh structure
    mesh = bpy.data.meshes.new("myMesh_mesh")
    mesh.from_pydata(vertsData, [], facesData)
    mesh.update()

    new_object = bpy.data.objects.new("myMesh_object", mesh)
    new_object.data = mesh

    scene = bpy.context.scene
    scene.objects.link(new_object)
    scene.objects.active = new_object
    new_object.select = True
    '''

    '''
    verts = []
    for i in range(0, len(tm_verts), 3):
        verts.append( (tm_verts[i], tm_verts[i+1], tm_verts[i+2]) )

    #mesh.edges.add(nbEdges)
    #print("# edges", nbEdges)
    #for i in range(nbEdges):
    #    mesh.edges[i].vertices = [edgesData[i*3],edgesData[i*3+1]]

    mesh.faces.add(nbFaces)
    print("# faces", nbFaces)
    for i in range(len(facesData)):
        mesh.faces[i].vertices = facesData[i]


    mesh = bpy.data.meshes.new("mesh")  # add a new mesh
    obj = bpy.data.objects.new("result_A", mesh)  # add a new object using the mesh

    scene = bpy.context.scene
    scene.objects.link(obj)  # put the object into the scene (link)
    scene.objects.active = obj  # set as the active object in the scene
    obj.select = True  # select object

    mesh = bpy.context.object.data
    bm = bmesh.new()

    for v in verts:
        bm.verts.new(v)  # add a new vert

    # make the bmesh the object's mesh
    bm.to_mesh(mesh)
    bm.free()  # always do this when finished
    '''

def register():
    bpy.utils.register_class(CSG)


def unregister():
    bpy.utils.unregister_class(CSG)


# This allows you to run the script directly from blenders text editor
# to test the addon without having to install it.
if __name__ == "__main__":
    register()
