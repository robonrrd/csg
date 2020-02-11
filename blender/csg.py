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

        # Blender objects
        clay_obj = None
        knife_obj = None
        for obj in scene.objects:
            if obj.type in ['MESH']:
                if clay_obj is None:
                    print("Using",obj.name,"as the clay")
                    clay_obj = obj
                elif knife_obj is None:
                    print("Using",obj.name,"as the knife")
                    knife_obj = obj
                else:
                    break

        clay = build_trimesh_from_obj(clay_obj)
        knife = build_trimesh_from_obj(knife_obj)

        if clay is not None and knife is not None:
            engine = pycsg.CSGEngine(clay, knife)
            A = pycsg.TriMesh()
            B = pycsg.TriMesh()
            engine.construct(pycsg.kDifference, True, A, B)

            build_blender_mesh_from_trimesh(A)
            build_blender_mesh_from_trimesh(B)

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


def build_blender_mesh_from_trimesh(trimesh):
    tm_verts = trimesh.vertices()
    print( type(tm_verts) )
    print( dir(tm_verts) )
    print( tm_verts.size() )
    #print("# verts", len(tm_verts))
    #print("  verts", tm_verts)

    #tm_faces = trimesh.faces()
    #print("# faces", len(tm_faces))
    #print("  faces", tm_faces)


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
