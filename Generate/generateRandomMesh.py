import bpy
import mathutils
import random
import math

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

import sys
argv = sys.argv
if "--" in argv:
    argv = argv[argv.index("--")+1:]
    if len(argv) > 0:
        if is_number(argv[0]):
            print("Using seed " + str(int(argv[0])))
            random.seed( int(argv[0]) )
    if len(argv) > 1:
        outputFolder = argv[1]+'/'



# Go to object mode if possible:
if bpy.context.object != None:
    if bpy.context.object.mode != 'OBJECT':
        bpy.ops.object.mode_set(mode='OBJECT')




# Clear Scene:
if "RandomMesh" in bpy.data.objects:
    bpy.data.objects["RandomMesh"].select_set(True)



if "Curve" in bpy.data.objects:
    bpy.data.objects["Curve"].select_set(True)




#for item in bpy.data.objects:
#    item.select = True
bpy.ops.object.delete()

radius = random.uniform(0.03,0.1)

vertices = []
vertices.append((0,radius,radius))
vertices.append((0,-radius,radius))
vertices.append((0,radius,-radius))
vertices.append((0,-radius,-radius))
vertices.append((0.05,radius,radius))
vertices.append((0.05,-radius,radius))
vertices.append((0.05,radius,-radius))
vertices.append((0.05,-radius,-radius))

edges = []
faces = []

faces.append((0,2,3,1))

faces.append((0,1,5,4))
faces.append((0,4,6,2))
faces.append((1,3,7,5))
faces.append((7,3,2,6))

faces.append((4,5,7,6))

mesh = bpy.data.meshes.new("RandomMesh")
mesh.from_pydata( vertices, edges, faces )
object = bpy.data.objects.new("RandomMesh", mesh )
scene = bpy.context.scene


#scene.objects.link(object) # was this in the older version of Blender
bpy.context.collection.objects.link(object)

bpy.ops.object.select_all(action='DESELECT')
# object.select = True
object.select_set(True)
#bpy.context.scene.objects.active = object
bpy.context.view_layer.objects.active = object

bpy.ops.object.mode_set(mode='EDIT')

#bpy.ops.object.editmode_toggle()
bpy.ops.mesh.select_all(action = 'DESELECT')


bpy.ops.object.mode_set(mode = 'OBJECT')
object.data.polygons[5].select = True
bpy.ops.object.mode_set(mode = 'EDIT')

extrusions = random.randint(1,4)

for i in range( 0, extrusions ):
    #bpy.ops.mesh.extrude_region_move(MESH_OT_extrude_region={"mirror":False}, TRANSFORM_OT_translate={"value":(0, 0, 0.05), "constraint_axis":(False, False, True), "constraint_orientation":'NORMAL', "mirror":False, "proportional":'DISABLED', "proportional_edit_falloff":'SMOOTH', "proportional_size":1, "snap":False, "snap_target":'CLOSEST', "snap_point":(0, 0, 0), "snap_align":False, "snap_normal":(0, 0, 0), "gpencil_strokes":False, "texture_space":False, "remove_on_cancel":False, "release_confirm":False})
    bpy.ops.mesh.extrude_region_move(TRANSFORM_OT_translate={"value":(0.05,0,0)})
    scale = random.uniform( 0.8,1.2 )
    bpy.ops.transform.resize(value=(scale,scale,scale));



bpy.ops.object.mode_set(mode = 'OBJECT')


for edge in mesh.edges:
    if random.randint(1,10) > 8:
        edge.crease = random.uniform(0,1)






#bpy.context.scene.cursor_location = (0.0, 0.0, 0.0)
bpy.context.scene.cursor.location = (0.0, 0.0, 0.0)
#bpy.ops.curve.primitive_bezier_curve_add(radius=1, view_align=False, enter_editmode=False)

# weight
w = 1


prevVector = mathutils.Vector((0,0,0))
pointList = [prevVector.copy()]
for i in range( 0, 10 ):
    prevVector += mathutils.Vector((0.05,random.uniform(-0.02,0.02), random.uniform(-0.02,0.02)))
    pointList.append( prevVector.copy() )



curvedata = bpy.data.curves.new(name="Curve", type='CURVE')
curvedata.dimensions = '3D'

curveObject = bpy.data.objects.new("Curve", curvedata)
curveObject.location = (0,0,0) #object origin


bpy.context.collection.objects.link(curveObject)
#bpy.context.scene.objects.link(curveObject)


polyline = curvedata.splines.new('NURBS')
polyline.points.add(len(pointList)-1)
for num in range(len(pointList)):
    polyline.points[num].co = pointList[num].to_4d()



polyline.order_u = len(polyline.points)-1
polyline.use_endpoint_u = True


modifierSubsurf = object.modifiers.new(name="Subsurf", type='SUBSURF')
#modifierSubsurf.levels = random.randint(3,4)
modifierSubsurf.levels = 3
modifierCurve = object.modifiers.new(name="Curve", type='CURVE')
modifierCurve.object = curveObject
bpy.ops.object.modifier_apply(modifier="Subsurf")
bpy.ops.object.modifier_apply(modifier="Curve")

# Randomly rotate the object:
eul = mathutils.Euler((random.uniform(0,2*math.pi),random.uniform(0,2*math.pi),random.uniform(0,2*math.pi
)), 'XYZ')
object.rotation_mode = "QUATERNION"
object.rotation_quaternion = eul.to_quaternion()
bpy.ops.object.transform_apply(location=False, rotation=True, scale=False)

# Move the object (and the curve) to fit into a predefined bounding box:
local_bbox_center = 0.125 * sum((mathutils.Vector(b) for b in object.bound_box), mathutils.Vector())
origOffset = -local_bbox_center     # Offset to move object into world origin
object.location = origOffset
curveObject.location = origOffset

# Offset to randomly place object inside bounding box:
maxOffset = (mathutils.Vector((0.3,0.3,0.3)) - object.dimensions)*0.5;
randOffset = mathutils.Vector(( random.uniform(-maxOffset.x,maxOffset.x),
                    random.uniform(-maxOffset.y,maxOffset.y),
                    random.uniform(-maxOffset.z,maxOffset.z) ))



print( maxOffset )
print( randOffset )
object.location += randOffset
curveObject.location += randOffset



filepath=outputFolder+str('randomMesh.stl')
# print("Exporting .stl to " + filepath);
bpy.ops.export_mesh.stl(filepath=filepath, check_existing=False, filter_glob="*.stl", ascii=True, use_mesh_modifiers=True, axis_forward='Y', axis_up='Z', global_scale=1.0)

# bpy.ops.wm.quit_blender()
