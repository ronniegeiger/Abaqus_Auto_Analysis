import os
import time
import random
from abaqus import *
from abaqusConstants import *
from Abaqus_Python_Functions import *
import numpy as np
import regionToolset
import __main__
from section import * 
import part
import material
import assembly
import step
import interaction
import load
import mesh
import job
import sketch
import visualization
import xyPlot
import connectorBehavior
import odbAccess

attemps = 40
set_size = 20

data = np.zeros((attemps,6))
data_save_path = r'D:\Desktop\Abaqus-python\CollectedData'

path = r'D:\Desktop\Abaqus-python\Abaqus\Data_41-80'

for k in range(attemps/set_size):
    i = k + 1
    Mdb()
    os.mkdir(r'{}\DataSet_{}'.format(path, i))
    for m in range(set_size):
        n = m + 1
        os.chdir(r'{}\DataSet_{}'.format(path, i))
        # Variables
        myRadius = random.randrange(30, 100)
        myLength = random.randrange(50, 100)
        myThickness = round(random.uniform(1.0, 2.0), 2)
        myOrientationAngle = 0.0
        BoundaryPerturbationAngle = round(random.uniform(0.0, 0.5), 2) * np.pi / 180.0
        BPangleindegree = round(BoundaryPerturbationAngle * 180.0 / np.pi, 3)
        
        LonElement = 100
        CirElement = 100
        
        na = myLength
        nc = 240
        
        myModelName = "Buckling_Analysis_{}".format(set_size*(i-1) + n)  
        myCylinderPart = "PART-1"
        myPlatePart1 = "UpperPlate"
        myPlatePart2 = "BottomPlate"
        myPlateLength = 200
        
        MonitoringSetName = "Top_RP"
        CylinderGeometrySet = "Cylinder_Set"
        
        # Material Properties
        myMaterial = "Woven"
        density = 1.45E-009
        e11 = 59500.0
        e22 = 55800.0
        e33 = 5900.0
        nu12 = 0.064
        nu13 = 0.064
        nu23 = 0.32
        g12 = 3650.0
        g13 = 3650.0
        g23 = 3650.0
        x1t = 540.0
        x1c = 440.0
        x2t = 440.0
        x2c = 380.0
        s12 = 40.0
        s13 = 150.0
        lte = 45.8
        lce = 39.95
        tte = 45.8
        tce = 39.95
        damage_stabilization = 0.0001
        
        # Step
        TimePeriod = 0.0001
        StepName = "Buckling"
        
        # Boundary Condition
        CompDistance = 2
        TopCondSet = "Top_Cond"
        BotCondSet = "Bot_Cond"
        
        # Function
        print('{} has been created'.format(myModelName))
        print('Length = {}'.format(myLength))
        print('Radius = {}'.format(myRadius))
        
        data[set_size*(i-1)+m,0] = myRadius
        data[set_size*(i-1)+m,1] = myLength
        data[set_size*(i-1)+m,2] = myThickness
        data[set_size*(i-1)+m,3] = BPangleindegree
        data[set_size*(i-1)+m,4] = Create_Random_Geometric_Imperfection(LonElement, CirElement, myThickness, myLength, myRadius, myModelName)
        Create_Node_File(myModelName, myLength, na, nc, myRadius)
        Create_Loading_Imperfection_Node_File(path, myLength, myRadius, BPangleindegree, myModelName, na, nc, i)
        Create_Inp_File(path, i, myModelName)
        Create_Element_File(myLength, nc, myModelName)
        Import_Inp_File(path, i, myModelName)
        Delete_Element(myModelName)
        Create_Part(myPlateLength, myPlatePart1, myPlatePart2, myModelName)
        Create_Cylinder_Geometry_Set(myModelName, CylinderGeometrySet, myLength, nc)
        Create_Material_Data_Input(myModelName, myMaterial, density, e11, e22, e33, nu12, nu13, nu23, g12, g13, g23, x1t, x1c, x2t, x2c, s12, s13, lte, lce, tte, tce, damage_stabilization)
        Assign_Material_Section(myModelName, CylinderGeometrySet, myOrientationAngle, myMaterial, myThickness)
        Create_Assembly(myModelName, myCylinderPart, myPlatePart1, myPlatePart2, myLength)
        Create_Interaction(myModelName)
        Create_Step(myModelName, TimePeriod, StepName)
        Create_Condition_Set(myModelName, TopCondSet, BotCondSet)
        Create_Condition(myModelName, BotCondSet, TopCondSet, TimePeriod, CompDistance, StepName)
        Create_Mesh(myModelName, myPlatePart2, myPlatePart1)
        Create_Job(myModelName)   
        Submit_Job(myModelName)
        Open_ODB_and_Write_NodeSet_data_to_text(myModelName,StepName,'U',TopCondSet.upper(),2)
        Open_ODB_and_Write_NodeSet_data_to_text(myModelName,StepName,'RF',TopCondSet.upper(),2)
        data[set_size*(i-1)+m,5] = Open_ODB_and_Write_Max_Value_of_NodeSet_data_to_text(myModelName,StepName,'RF','N',TopCondSet.upper(),2)
    mdb.saveAs(pathName=r"{}\DataSet_{}\CAE_Set{}".format(path, i, i))

np.save(r'{}\Data41-80'.format(data_save_path),data)
f = np.load(r'{}\Data41-80.npy'.format(data_save_path))



