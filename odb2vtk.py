'''
@ The script is developed by Qingbin Liu under the supervision of Jie Liu
E-mail to Qingbin Liu: liuqingb@mail2.sysu.edu.cn; liuqingb@pku.edu.cn
       to discuss technical details. 
--------------------------------------------------------------------------------
       Github: Echeban/odb2vtk. Modifications to Function1 only.
       Function 2: kept unchanged as in 2016 commit main branch
--------------------------------------------------------------------------------
The script is to convert the data from Abaqus output database format 
to vtk file format for parallel visualization
Python script performs the following three major steps: 
    1) reading Abaqus output file according to the architecture of ODB; 
    2) data decomposition for parallel visualization; 
    3) writing VTK for Paraview.
'''

#import necessary modules to handle Abaqus output database, files and string
from odbAccess import *
from textRepr import *
from string import *
from time import *

#Function 1
def ConvertOdb2Vtk(filename = 'C:\SIMULIA\User\odb2vtk.txt'):  #Modify the default value of filename here to specify the default configuration file
    
    starttime = time()
    # Check if filename points to existing file
    if not os.path.isfile(filename): 
        print 'Parameter file "%s" not found'%(filename)
        sys.exit(2) 
    
    #read odb2vtk file to get parameter setting 
    odb2vtk = open(filename,'rt')
    read = odb2vtk.read()
    input = read.split("'")
    #get odb file's path
    odb_path = input[1]
    #get odb file name
    odbname = input[3]
    #get the output files' path
    vtk_path = input[5]
    #get the mesh type
    mesh_type = int(input[7])
    mesh_conner = 0
    if (mesh_type == 12):
        mesh_conner = 8
        mesh_name = "Hexahedron"
    if (mesh_type == 10):
        mesh_conner = 4 
        mesh_name = "Tetra"
    if (mesh_type == 9):
        mesh_conner = 4 # EJB: number of corner nodes
        mesh_name = "Quad"
    if (mesh_conner == 0):
        print "Mesh type error or unidentified"
        os._exit(0)
    #get the quantity of pieces to partition
    piecenum = int(input[9])
    #get the frame
    input_frame = input[11].split("-")
    input_frame = range(int(input_frame[0]),int(input_frame[1])+1)
    #get the step
    input_step = input[13].split(",")
    #get the instance
    input_instance = input[15].split(",")
    #end reding and close odb2vtk file
    odb2vtk.close()
    #display the reading result of odb2vtk file
    print "odb2vtk reading finished, time elapsed: ", time()-starttime
    print "Basic Information:"
    print "Model:",odbname,"; Mesh type:",mesh_name,"; Number of blocks:",piecenum
    print "Convert frames: ",input_frame[0]," to ",input_frame[-1]
    print "Step & Instance : ",str(input_step),", ",str(input_instance)
    
    #open an ODB ( Abaqus output database )
    odb = openOdb(os.path.join(odb_path,odbname)+'.odb',readOnly=True)
    print "ODB opened"

    #access geometry and topology information ( odb->rootAssembly->instances->(nodes, elements) )
    rootassembly = odb.rootAssembly
    instance = rootassembly.instances
    #access attribute information
    step = odb.steps
    #get instance & step information : Quantity and all names
    allinstancestr = str(instance)
    autoins = allinstancestr.split("'")
    inslen = len(autoins)/4
    instance_N = range(0,inslen)
    allstepstr = str(step)
    autostep = allstepstr.split("'")
    steplen = len(autostep)/4
    step_N = range(0,steplen)
    
    for i in input_step:
        if(steplen < int(i)):
            print "input step exceeds the range of steps"
            os._exit(0)
    for i in input_instance:
        if(inslen < int(i)):
            print "input instance exceeds the range of instances"
            os._exit(0)
        

    #step cycle
    for step_i in input_step:
        n = int(step_i)*4+1
        stepname = autostep[n]
        print "Step: ",stepname
        #instance cycle
        for ins_i in input_instance:
            n = int(ins_i)*4+1
            instancename = autoins[n]
            print "Instance: ",instancename
            
            #access nodes & elements
            node = instance[instancename].nodes
            element = instance[instancename].elements
            n_nodes = len(node)
            n_elements = len(element)
            #access attribute(fieldOutputs) information
            frame = step[stepname].frames
            
            
            #compute the number of element of each block
            p_elements = n_elements/piecenum + 1
            lp_elements = n_elements - (p_elements*(piecenum-1))  #last block
            
            
            #match nodes' label and its order in sequence (for empty nodes in tetra mesh)
            MLN = node[n_nodes-1].label
            TOTAL=[]
            #read node in sequence, and get the largest label of node(non-empty) 
            #MLN is the max label of nodeset
            for i in node:
                TOTAL.append(i.label)
                if(i.label > MLN):
                    MLN = i.label
            #match (the key)
            L=[]
            n = 0
            for i in range(MLN): 
                L.append(0)
            for i in TOTAL:
                L[i-1] = n
                n += 1

            
            #frame cycle
            for i_frame in input_frame:
                
                #Detect whether the input frame is out of range
                try:
                    TRY = odb.steps[stepname].frames[int(i_frame)]
                except:
                    print "* warning * input frame exceeds the range of frames" 
                    # os._exit(0) #just break the loop, no need to kill the Python shell
                    break
                
                #Access a frame
                N_Frame = odb.steps[stepname].frames[int(i_frame)]
                print "Frame:",i_frame
                
                #create array for store result data temporarily
                # Vector-U,A,V,RF 
                L0=[] 
                # Tensors-S
                L1=[]
                # Tensors-LE,E
                L2=[]
                L5=[] #strain E
                # Tensors-PE
                L3=[]
                # Scalars-PEEQ
                L4=[]
                for i in range(MLN): 
                    L0.append([0,0,0,0,0,0,0,0,0,0,0,0])
                    L1.append([1,0,0,0,0,0,0,0,0,0,0,0,0,0])
                    L2.append([1,0,0,0,0,0,0,0,0]) 
                    L3.append([1,0,0,0,0,0,0,0,0])
                    L4.append([1,0]) #1, prevents /zero when no data for this keyword
                    L5.append([1,0,0,0,0,0,0,0,0])
                
                print "Reading U, A, V, RF ......"
                time1 = time()
                #Access Spatial displacement
                try:
                    displacement = N_Frame.fieldOutputs['U']
                    fieldValues = displacement.values
                    for valueX in fieldValues :
                        i = valueX.nodeLabel
                        L0[i-1][0] = round(valueX.data[0],37)#round() : 1E-37 < ParaView-values < 1E37
                        L0[i-1][1] = round(valueX.data[1],37)
                        L0[i-1][2] = round(valueX.data[2],37)
                    UDataExists=True
                    
                except KeyError:
                    print "No data for keyword : U"
                    UdataExists=False
                
                # Access Spatial acceleration
                try:
                    acceleration = N_Frame.fieldOutputs['A']
                    fieldValues = acceleration.values
                    for valueX in fieldValues :
                        i = valueX.nodeLabel
                        L0[i-1][3] = valueX.data[0]
                        L0[i-1][4] = valueX.data[1]
                        L0[i-1][5] = valueX.data[2]
                    ADataExists=True
                    
                except KeyError:
                    print "No data for keyword : A"
                    ADataExists=False
                
                # Access Spatial velocity
                try:
                    velocity = N_Frame.fieldOutputs['V']
                    fieldValues = velocity.values
                    for valueX in fieldValues :
                        i = valueX.nodeLabel
                        L0[i-1][6] = valueX.data[0]
                        L0[i-1][7] = valueX.data[1]
                        L0[i-1][8] = valueX.data[2]
                    VDataExists=True
                    
                except KeyError:
                    print "No data for keyword : V"
                    VDataExists=False
                
                #Access Reaction force
                try:
                    Reaction_force = N_Frame.fieldOutputs['RF']
                    fieldValues = Reaction_force.values 
                    for valueX in fieldValues :
                        i = valueX.nodeLabel
                        L0[i-1][9] = valueX.data[0]
                        L0[i-1][10] = valueX.data[1]
                        L0[i-1][11] = valueX.data[2]    
                    RFDataExists=True
                    print "Time elapsed: ", time() - time1, "s"
                except KeyError:
                    print "No data for keyword : RF"
                    RFDataExists=False
                
                print "Reading Stress ......"
                time1 = time()
                #Access Stress components
                try:
                    Stress = N_Frame.fieldOutputs['S']
                    node_Stress = Stress.getSubset(position=ELEMENT_NODAL)
                    fieldValues = node_Stress.values
                    for valueX in fieldValues :
                        L1[valueX.nodeLabel-1][0] += 1
                        L1[valueX.nodeLabel-1][1] += valueX.data[0]
                        L1[valueX.nodeLabel-1][2] += valueX.data[1]
                        L1[valueX.nodeLabel-1][3] += valueX.data[2]
                        L1[valueX.nodeLabel-1][4] += valueX.data[3]
                        if mesh_type==10 or mesh_type==12: #only continuum elements have 6 stress values
                            L1[valueX.nodeLabel-1][5] += valueX.data[4]
                            L1[valueX.nodeLabel-1][6] += valueX.data[5]
                        
                        L1[valueX.nodeLabel-1][7] += valueX.mises
                        L1[valueX.nodeLabel-1][8] += valueX.maxPrincipal
                        L1[valueX.nodeLabel-1][9] += valueX.midPrincipal
                        L1[valueX.nodeLabel-1][10] += valueX.minPrincipal
                        L1[valueX.nodeLabel-1][11] += valueX.press
                        L1[valueX.nodeLabel-1][12] += valueX.tresca
                        L1[valueX.nodeLabel-1][13] += valueX.inv3
                    SDataExists=True
                    
                except KeyError:
                    print "No data for keyword : S"
                    SDataExists=False
                
                # can first ave
                print "Time elapsed: ", time() - time1, "s"
                print "Reading Logarithmic strain ......"
                time1 = time()
                # Logarithmic strain components
                try:
                    Logarithmic_strain = N_Frame.fieldOutputs['LE']
                    node_Logarithmic_strain = Logarithmic_strain.getSubset(position=ELEMENT_NODAL)
                    fieldValues = node_Logarithmic_strain.values
                    for valueX in fieldValues :
                        L2[valueX.nodeLabel-1][0] += 1
                        L2[valueX.nodeLabel-1][1] += valueX.data[0]
                        L2[valueX.nodeLabel-1][2] += valueX.data[1]
                        L2[valueX.nodeLabel-1][3] += valueX.data[2]
                        L2[valueX.nodeLabel-1][4] += valueX.data[3]
                        if mesh_type==10 or mesh_type==12: #only cont elems have 6 stress values
                            L2[valueX.nodeLabel-1][5] += valueX.data[4]
                            L2[valueX.nodeLabel-1][6] += valueX.data[5]
                        
                        L2[valueX.nodeLabel-1][7] += valueX.maxPrincipal
                        L2[valueX.nodeLabel-1][8] += valueX.minPrincipal
                    LEDataExists=True
                    print "Time elapsed: ", time() - time1, "s"
                except KeyError:
                    print "No data for keyword : LE"
                    LEDataExists=False
                
                print "Reading Plastic strain ......"
                time1 = time()
                #Plastic strain components
                try:
                    Plastic_strain = N_Frame.fieldOutputs['PE']
                    node_Plastic_strain = Plastic_strain.getSubset(position=ELEMENT_NODAL)
                    fieldValues = node_Plastic_strain.values    
                    for valueX in fieldValues :
                        L3[valueX.nodeLabel-1][0] += 1
                        L3[valueX.nodeLabel-1][1] += valueX.data[0]
                        L3[valueX.nodeLabel-1][2] += valueX.data[1]
                        L3[valueX.nodeLabel-1][3] += valueX.data[2]
                        L3[valueX.nodeLabel-1][4] += valueX.data[3]
                        if mesh_type==10 or mesh_type==12: #only cont elems have 6 stress values
                            L3[valueX.nodeLabel-1][5] += valueX.data[4]
                            L3[valueX.nodeLabel-1][6] += valueX.data[5]
                        
                        L3[valueX.nodeLabel-1][7] += valueX.maxPrincipal
                        L3[valueX.nodeLabel-1][8] += valueX.minPrincipal
                    PEDataExists=True
                    print "Time elapsed: ", time() - time1, "s"
                except KeyError:
                    print "No data for keyword : PE"
                    PEDataExists=False
                
                print "Reading Equivalent plastic strain ......"
                time1 = time()
                #Equivalent plastic strain
                try:
                    Equivalent_plastic_strain = N_Frame.fieldOutputs['PEEQ']
                    node_Equivalent_plastic_strain = Equivalent_plastic_strain.getSubset(position=ELEMENT_NODAL)
                    fieldValues = node_Equivalent_plastic_strain.values
                    for valueX in fieldValues :
                        L4[valueX.nodeLabel-1][0] += 1
                        L4[valueX.nodeLabel-1][1] += valueX.data
                    PEEQDataExists=True
                    print "Time elapsed: ", time() - time1, "s" 
                except KeyError:
                    print "No data for keyword : PEEQ"
                    PEEQDataExists=False
                
                print "Reading total strain ......"
                time1 = time()
                # strain components
                try:
                    Strain = N_Frame.fieldOutputs['E']
                    node_strain = Strain.getSubset(position=ELEMENT_NODAL)
                    fieldValues = node_strain.values
                    for valueX in fieldValues :
                        L5[valueX.nodeLabel-1][0] += 1
                        L5[valueX.nodeLabel-1][1] += valueX.data[0]
                        L5[valueX.nodeLabel-1][2] += valueX.data[1]
                        L5[valueX.nodeLabel-1][3] += valueX.data[2]
                        L5[valueX.nodeLabel-1][4] += valueX.data[3]
                        if mesh_type==10 or mesh_type==12: #only cont elems have 6 stress values
                            L5[valueX.nodeLabel-1][5] += valueX.data[4]
                            L5[valueX.nodeLabel-1][6] += valueX.data[5]
                        
                        L5[valueX.nodeLabel-1][7] += valueX.maxPrincipal
                        L5[valueX.nodeLabel-1][8] += valueX.minPrincipal
                    EDataExists=True
                    print "Time elapsed: ", time() - time1, "s"
                except KeyError:
                    print "No data for keyword : E"
                    EDataExists=False
                
                '''============================================================'''
                
                print "Partitioning model and writing vtk files ......"
                #piece cycle, to partition the model and create each piece for vtk files      
                for pn in range(piecenum):
                    time1 = time()
                    print "frame:",i_frame,"; block:",pn
                    #Reorganization
                    #Control&Storage
                    #estimate whether the node has already existed
                    stg_p = []
                    #store the reorganized node for element
                    stg_e = []
                    #store the reorganized node for node
                    stg_n = []
                    for i in range(MLN):
                        stg_p.append(-1)
                    nodecount = 0
                    #reorganize the node and element (reconstruct the mesh)
                    if(pn == piecenum-1):
                        M = range(pn*p_elements,n_elements)
                    else:
                        M = range(pn*p_elements,(pn+1)*p_elements)
                    for i in M:
                        for j in range(mesh_conner):
                            k = element[i].connectivity[j] - 1
                            if(stg_p[k] < 0): 
                                stg_p[k] = nodecount
                                stg_n.append(L[k]) 
                                stg_e.append(nodecount)
                                nodecount += 1
                            else:
                                stg_e.append(stg_p[k])
                    #compute point quantity
                    n_reop = len(stg_n)
                    reop_N = range(0,len(stg_n))


                    #create and open a VTK(.vtu) files
                    if(piecenum > 1):
                        outfile = open (os.path.join(vtk_path,odbname)+'_'+stepname+'_'+instancename+'f%03d'%int(i_frame)+' '+'p'+str(pn)+'.vtu','w')
                    if(piecenum == 1):
                        outfile = open (os.path.join(vtk_path,odbname)+'_'+stepname+'_'+instancename+'f%03d'%int(i_frame)+'.vtu','w')
                    
                    #<VTKFile>, including the type of mesh, version, and byte_order
                    outfile.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'+'\n')
                    #<UnstructuredGrid>
                    outfile.write('<UnstructuredGrid>'+'\n')
                    #<Piece>, including the number of points and cells
                    if(pn == piecenum-1):
                        outfile.write('<Piece NumberOfPoints="'+str(n_reop)+'"'+' '+'NumberOfCells="'+str(lp_elements)+'">'+'\n')
                    else:
                        outfile.write('<Piece NumberOfPoints="'+str(n_reop)+'"'+' '+'NumberOfCells="'+str(p_elements)+'">'+'\n')

                    
                    print "Writing Nodes ......"
                    #<Points> Write nodes into vtk files
                    displacement = N_Frame.fieldOutputs['U']
                    fieldValues = displacement.values
                    outfile.write('<Points>'+'\n')
                    outfile.write('<DataArray type="Float64" NumberOfComponents="3" format="ascii">'+'\n')
                    for i in reop_N:
                        nt = stg_n[i]
                        k = node[stg_n[i]].label-1
                        X,Y,Z = node[nt].coordinates[0]+L0[k][0],node[nt].coordinates[1]+L0[k][1],node[nt].coordinates[2]+L0[k][2]
                        outfile.write(' '+'%11.8e'%X+'  '+'%11.8e'%Y+'  '+'%11.8e'%Z+'\n')          
                    outfile.write('</DataArray>'+'\n')
                    outfile.write('</Points>'+'\n')
                    #</Points>


                    print "Writing Results data ......"
                    #<PointData> Write results data into vtk files
                    outfile.write("<"+"PointData"+" "+"Tensors="+'"'+"Stress_Components,Plastic_strain_components"+'"'\
                    +" "+"Vectors="+'"'+"Spatial_displacement,Reaction_force"+'"'\
                    +" "+"Scalars="+'"'+"Equivalent_plastic_strain,Stress_Mises,Stress_Max_Principal,Stress_Mid_Principal,Stress_Min_Principal,Stress_Pressure,Stress_Tresca,Stress_Third_Invariant,Plastic_strain_Max_Principal,Plastic_strain_Min_Principal"+'"'+">"+'\n')
                    
                    #Stress components S, <DataArray>
                    outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Components"+'"'+" "+"NumberOfComponents="+'"'+"9"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
                    for i in reop_N:
                        k = node[stg_n[i]].label-1
                        XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ = L1[k][1]/L1[k][0],L1[k][4]/L1[k][0],L1[k][6]/L1[k][0],L1[k][4]/L1[k][0],L1[k][2]/L1[k][0],L1[k][5]/L1[k][0],L1[k][6]/L1[k][0],L1[k][5]/L1[k][0],L1[k][3]/L1[k][0]
                        outfile.write('%11.8e'%XX+' '+'%11.8e'%XY+' '+'%11.8e'%XZ+' '+'%11.8e'%YX+' '+'%11.8e'%YY+' '+'%11.8e'%YZ+' '+'%11.8e'%ZX+' '+'%11.8e'%ZY+' '+'%11.8e'%ZZ+'\n')
                    outfile.write("</DataArray>"+'\n')
                    #</DataArray>

                    #Strain components E, <DataArray>
                    outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Strain_components"+'"'+" "+"NumberOfComponents="+'"'+"9"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
                    for i in reop_N:
                        k = node[stg_n[i]].label-1
                        XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ = L5[k][1]/L5[k][0],L5[k][4]/L5[k][0],L5[k][6]/L5[k][0],L5[k][4]/L5[k][0],L5[k][2]/L5[k][0],L5[k][5]/L5[k][0],L5[k][6]/L5[k][0],L5[k][5]/L5[k][0],L5[k][3]/L5[k][0]
                        outfile.write('%11.8e'%XX+' '+'%11.8e'%XY+' '+'%11.8e'%XZ+' '+'%11.8e'%YX+' '+'%11.8e'%YY+' '+'%11.8e'%YZ+' '+'%11.8e'%ZX+' '+'%11.8e'%ZY+' '+'%11.8e'%ZZ+'\n')
                    outfile.write("</DataArray>"+'\n')
                    #</DataArray>
                    
                    #Logarithmic strain components, <DataArray>
                    outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Logarithmic_strain_components"+'"'+" "+"NumberOfComponents="+'"'+"9"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
                    for i in reop_N:
                        k = node[stg_n[i]].label-1
                        XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ = L2[k][1]/L2[k][0],L2[k][4]/L2[k][0],L2[k][6]/L2[k][0],L2[k][4]/L2[k][0],L2[k][2]/L2[k][0],L2[k][5]/L2[k][0],L2[k][6]/L2[k][0],L2[k][5]/L2[k][0],L2[k][3]/L2[k][0]
                        outfile.write('%11.8e'%XX+' '+'%11.8e'%XY+' '+'%11.8e'%XZ+' '+'%11.8e'%YX+' '+'%11.8e'%YY+' '+'%11.8e'%YZ+' '+'%11.8e'%ZX+' '+'%11.8e'%ZY+' '+'%11.8e'%ZZ+'\n')
                    outfile.write("</DataArray>"+'\n')
                    #</DataArray>
                    
                    #Plastic strain components, <DataArray>
                    outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Plastic_strain_components"+'"'+" "+"NumberOfComponents="+'"'+"9"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
                    for i in reop_N:
                        k = node[stg_n[i]].label-1
                        XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ = L3[k][1]/L3[k][0],L3[k][4]/L3[k][0],L3[k][6]/L3[k][0],L3[k][4]/L3[k][0],L3[k][2]/L3[k][0],L3[k][5]/L3[k][0],L3[k][6]/L3[k][0],L3[k][5]/L3[k][0],L3[k][3]/L3[k][0]
                        outfile.write('%11.8e'%XX+' '+'%11.8e'%XY+' '+'%11.8e'%XZ+' '+'%11.8e'%YX+' '+'%11.8e'%YY+' '+'%11.8e'%YZ+' '+'%11.8e'%ZX+' '+'%11.8e'%ZY+' '+'%11.8e'%ZZ+'\n')
                    outfile.write("</DataArray>"+'\n')
                    #</DataArray>

                    #Spatial displacement, <DataArray>
                    outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Spatial_displacement"+'"'+" "+"NumberOfComponents="+'"'+"3"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
                    for i in reop_N:
                        k = node[stg_n[i]].label-1
                        X,Y,Z = L0[k][0],L0[k][1],L0[k][2]
                        outfile.write('%11.8e'%X+' '+'%11.8e'%Y+' '+'%11.8e'%Z+'\n')
                    outfile.write("</DataArray>"+'\n')
                    #</DataArray>
        
                    #Spatial acceleration, <DataArray>
                    outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Spatial_acceleration"+'"'+" "+"NumberOfComponents="+'"'+"3"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
                    for i in reop_N:
                        k = node[stg_n[i]].label-1
                        X,Y,Z = L0[k][3],L0[k][4],L0[k][5]
                        outfile.write('%11.8e'%X+' '+'%11.8e'%Y+' '+'%11.8e'%Z+'\n')
                    outfile.write("</DataArray>"+'\n')
                    #</DataArray>

                    #Spatial velocity, <DataArray>
                    outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Spatial_velocity"+'"'+" "+"NumberOfComponents="+'"'+"3"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
                    for i in reop_N:
                        k = node[stg_n[i]].label-1
                        X,Y,Z = L0[k][6],L0[k][7],L0[k][8]
                        outfile.write('%11.8e'%X+' '+'%11.8e'%Y+' '+'%11.8e'%Z+'\n')
                    outfile.write("</DataArray>"+'\n')  
                    #</DataArray>
        
                    #Reaction force
                    outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Reaction_force"+'"'+" "+"NumberOfComponents="+'"'+"3"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
                    for i in reop_N:
                        k = node[stg_n[i]].label-1
                        X,Y,Z = L0[k][9],L0[k][10],L0[k][11]
                        outfile.write('%11.8e'%X+' '+'%11.8e'%Y+' '+'%11.8e'%Z+'\n')
                    outfile.write("</DataArray>"+'\n')  
                    #</DataArray>
                    
                    #Equivalent plastic strain, <DataArray>
                    outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Equivalent_plastic_strain"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
                    for i in reop_N:
                        k = node[stg_n[i]].label-1
                        X = L4[k][1]/L4[k][0]
                        outfile.write('%11.8e'%X+'\n')
                    outfile.write('</DataArray>'+'\n')
                    #</DataArray>

                    #Stress Mises, <DataArray>
                    outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Mises"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
                    for i in reop_N:
                        k = node[stg_n[i]].label-1
                        X = L1[k][7]/L1[k][0]
                        outfile.write('%11.8e'%X+'\n')
                    outfile.write('</DataArray>'+'\n')
                    #</DataArray>

                    #Stress Max.Principal, <DataArray>
                    outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Max_Principal"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
                    for i in reop_N:
                        k = node[stg_n[i]].label-1
                        X = L1[k][8]/L1[k][0]
                        outfile.write('%11.8e'%X+'\n')
                    outfile.write('</DataArray>'+'\n')  
                    #</DataArray>
        
                    #Stress Mid.Principal, <DataArray>
                    outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Mid_Principal"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
                    for i in reop_N:
                        k = node[stg_n[i]].label-1
                        X = L1[k][9]/L1[k][0]
                        outfile.write('%11.8e'%X+'\n')
                    outfile.write('</DataArray>'+'\n')  
                    #</DataArray>
            
                    #Stress Min.Principal, <DataArray>
                    outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Min_Principal"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
                    for i in reop_N:
                        k = node[stg_n[i]].label-1
                        X = L1[k][10]/L1[k][0]
                        outfile.write('%11.8e'%X+'\n')
                    outfile.write('</DataArray>'+'\n')  
                    #</DataArray>
            
                    #Stress Pressure, <DataArray>
                    outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Pressure"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
                    for i in reop_N:
                        k = node[stg_n[i]].label-1
                        X = L1[k][11]/L1[k][0]
                        outfile.write('%11.8e'%X+'\n')
                    outfile.write('</DataArray>'+'\n')  
                    #</DataArray>

                    #Stress Tresca, <DataArray>
                    outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Tresca"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
                    for i in reop_N:
                        k = node[stg_n[i]].label-1
                        X = L1[k][12]/L1[k][0]
                        outfile.write('%11.8e'%X+'\n')
                    outfile.write('</DataArray>'+'\n')  
                    #</DataArray>
            
                    #Stress Third_Invariant, <DataArray>
                    outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Third_Invariant"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
                    for i in reop_N:
                        k = node[stg_n[i]].label-1
                        X = L1[k][13]/L1[k][0]
                        outfile.write('%11.8e'%X+'\n')
                    outfile.write('</DataArray>'+'\n')  
                    #</DataArray>
        
                    #Logarithmic_strain_Max_Principal, <DataArray>
                    outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Logarithmic_strain_Max_Principal"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
                    for i in reop_N:
                        k = node[stg_n[i]].label-1
                        X = L2[k][7]/L2[k][0]
                        outfile.write('%11.8e'%X+'\n')
                    outfile.write('</DataArray>'+'\n')  
                    #</DataArray>
        
                    #Logarithmic strain Min.Principal, <DataArray>
                    outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Logarithmic_strain_Min_Principal"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
                    for i in reop_N:
                        k = node[stg_n[i]].label-1
                        X = L2[k][8]/L2[k][0]
                        outfile.write('%11.8e'%X+'\n')
                    outfile.write('</DataArray>'+'\n')  
                    #</DataArray>'''
        
                    #Plastic strain Max.Principal, <DataArray>
                    outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Plastic_strain_Max_Principal"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
                    for i in reop_N:
                        k = node[stg_n[i]].label-1
                        X = L3[k][7]/L3[k][0]
                        outfile.write('%11.8e'%X+'\n')
                    outfile.write('</DataArray>'+'\n')  
                    #</DataArray>
        
                    #Plastic strain Min.Principal, <DataArray>
                    outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Plastic_strain_Min_Principal"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
                    for i in reop_N:
                        k = node[stg_n[i]].label-1
                        X = L3[k][8]/L3[k][0]
                        outfile.write('%11.8e'%X+'\n')
                    outfile.write('</DataArray>'+'\n')  
                    #</DataArray>
                    
                    outfile.write("</PointData>"+'\n')
                    #</PointData>
                    
                    
                    print "Writing Cells ......"
                    #<Cells> Write cells into vtk files
                    outfile.write('<Cells>'+'\n')
                    #Connectivity
                    outfile.write('<DataArray type="Int32" Name="connectivity" format="ascii">'+'\n')
                    if (mesh_type == 12):
                        for i in range(len(stg_e)/8):
                            outfile.write(str(stg_e[i*8])+' '+str(stg_e[i*8+1])+' '+str(stg_e[i*8+2])+' '+str(stg_e[i*8+3])+' '+str(stg_e[i*8+4])+' '+str(stg_e[i*8+5])+' '+str(stg_e[i*8+6])+' '+str(stg_e[i*8+7])+'\n')
                    if (mesh_type == 10):
                        for i in range(len(stg_e)/4):
                            outfile.write(str(stg_e[i*4])+' '+str(stg_e[i*4+1])+' '+str(stg_e[i*4+2])+' '+str(stg_e[i*4+3])+'\n')
                    if (mesh_type == 9):
                        for i in range(len(stg_e)/4):
                            outfile.write(str(stg_e[i*4])+' '+str(stg_e[i*4+1])+' '+str(stg_e[i*4+2])+' '+str(stg_e[i*4+3])+'\n') #EJB
                    outfile.write('</DataArray>'+'\n')
                    #Offsets
                    outfile.write('<DataArray type="Int32" Name="offsets" format="ascii">'+'\n')
                    for i in range(len(stg_e)/mesh_conner):
                        outfile.write(str(i*mesh_conner+mesh_conner)+'\n')
                    outfile.write('</DataArray>'+'\n')
                    #Type
                    outfile.write('<DataArray type="UInt8" Name="types" format="ascii">'+'\n')
                    for i in range(len(stg_e)/mesh_conner):
                        outfile.write(str(mesh_type)+'\n')
                    outfile.write('</DataArray>'+'\n')
                    outfile.write('</Cells>'+'\n')
                    #</Cells>
        
        
                    #</Piece>
                    outfile.write('</Piece>'+'\n')
                    #</UnstructuredGrid>
                    outfile.write('</UnstructuredGrid>'+'\n')
                    #</VTKFile>
                    outfile.write('</VTKFile>'+'\n')
                
                    outfile.close()
                    print "Time elapsed: ", time() - time1, "s" 
                
                '''====================================================================='''
                print "Creating .pvtu file for frame ", i_frame," ......"
                #create .pvtu files for parallel visualization
                if ( piecenum > 1 ):
                    outfile = open (os.path.join(vtk_path,odbname)+'_'+stepname+'_'+'f%03d'%int(i_frame)+'.pvtu','w')
                    
                    #write the basic information for .pvtu files
                    outfile.write('<?xml version="1.0"?>'+'\n')
                    outfile.write('<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">'+'\n')
                    outfile.write("<PUnstructuredGrid GhostLevel="+'"'+str(piecenum)+'"'+">"+'\n')
                    #pointdata
                    outfile.write("<"+"PPointData"+" "+"Tensors="+'"'+"Stress_Components,Logarithmic_strain_components,Plastic_strain_components"+'"'\
                    +" "+"Vectors="+'"'+"Spatial_displacement,Spatial_acceleration,Spatial_velocity,Reaction_force"+'"'\
                    +" "+"Scalars="+'"'+"Equivalent_plastic_strain,Stress_Mises,Stress_Max_Principal,Stress_Mid_Principal,Stress_Min_Principal,Stress_Pressure,Stress_Tresca,Stress_Third_Invariant,Logarithmic_strain_Max_Principal,Logarithmic_strain_Min_Principal,Plastic_strain_Max_Principal,Plastic_strain_Min_Principal"+'"'+">"+'\n')
                    
                    outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Components"+'"'+" "+"NumberOfComponents="+'"'+"9"+'"'+" "+"/>"+'\n')
                    outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Logarithmic_strain_components"+'"'+" "+"NumberOfComponents="+'"'+"9"+'"'+" "+"/>"+'\n')
                    outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Plastic_strain_components"+'"'+" "+"NumberOfComponents="+'"'+"9"+'"'+" "+"/>"+'\n')
                    outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Spatial_displacement"+'"'+" "+"NumberOfComponents="+'"'+"3"+'"'+" "+"/>"+'\n')
                    outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Spatial_acceleration"+'"'+" "+"NumberOfComponents="+'"'+"3"+'"'+" "+"/>"+'\n')
                    outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Spatial_velocity"+'"'+" "+"NumberOfComponents="+'"'+"3"+'"'+" "+"/>"+'\n')
                    outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Reaction_force"+'"'+" "+"NumberOfComponents="+'"'+"3"+'"'+" "+"/>"+'\n')
                    outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Equivalent_plastic_strain"+'"'+" "+"/>"+'\n')
                    outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Mises"+'"'+" "+"/>"+'\n')
                    outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Max_Principal"+'"'+" "+"/>"+'\n')
                    outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Mid_Principal"+'"'+" "+"/>"+'\n')
                    outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Min_Principal"+'"'+" "+"/>"+'\n')
                    outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Pressure"+'"'+" "+"/>"+'\n')
                    outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Tresca"+'"'+" "+"/>"+'\n')
                    outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Third_Invariant"+'"'+" "+"/>"+'\n')
                    outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Logarithmic_strain_Max_Principal"+'"'+" "+"/>"+'\n')
                    outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Logarithmic_strain_Min_Principal"+'"'+" "+"/>"+'\n')
                    outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Plastic_strain_Max_Principal"+'"'+" "+"/>"+'\n')
                    outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Plastic_strain_Min_Principal"+'"'+" "+"/>"+'\n')
                    outfile.write("</PPointData>"+'\n')
                    
                    #points
                    outfile.write("<PPoints>"+'\n')
                    outfile.write("<PDataArray type="+'"'+"Float64"+'"'+" "+"NumberOfComponents="+'"'+"3"+'"'+"/>"+'\n')
                    outfile.write("</PPoints>"+'\n')
                    
                    #write the path of each piece for reading it through the .pvtu file 
                    for pn in range(piecenum):
                        outfile.write("<Piece Source="+'"'+odbname+'_'+stepname+'_'+instancename+'f%03d'%int(i_frame)+' '+'p'+str(pn)+'.vtu'+'"'+"/>"+'\n')
                    
                    outfile.write("</PUnstructuredGrid>"+'\n')
                    outfile.write("</VTKFile>")

                    outfile.close() 
            
            odb.close()

    print "Total time elapsed: ", time() - starttime, "s"

#Function 2	
def ConvertOdb2VtkP(Odbpath = ' ',Odbname = ' ',vtkpath = ' ',MeshType = ' ',Piecenum = ' ', BeginFrame = ' ', EndFrame =' ', Steps =' ', Instances = ' '):
	
	starttime = time()
	
	#get odb file's path
	odb_path = Odbpath
	#get odb file name
	odbname = Odbname
	#get the output files' path
	vtk_path = vtkpath
	#get the mesh type
	mesh_type = int(MeshType)
	mesh_conner = 0
	if (mesh_type == 12):
		mesh_conner = 8
		mesh_name = "Hexahedron"
	if (mesh_type == 10):
		mesh_conner = 4
		mesh_name = "Tetra"
	if (mesh_conner == 0):
		print "Mesh type error or unidentified"
		os._exit(0)
	#get the quantity of pieces to partition
	piecenum = int(Piecenum)
	#get the frame
	input_frame = range(int(BeginFrame),int(EndFrame)+1)
	#get the step
	input_step = Steps.split(",")
	#get the instance
	input_instance = Instances.split(",")

	#display the read parameter
	print "Basic Information:"
	print "Model:",odbname,"; Mesh type:",mesh_name,"; Number of blocks:",piecenum
	print "Convert frames: ",input_frame[0]," to ",input_frame[-1]
	print "Step & Instance : ",str(input_step),", ",str(input_instance)
	
	#open an ODB ( Abaqus output database )
	odb = openOdb(os.path.join(odb_path,odbname)+'.odb',readOnly=True)
	print "ODB opened"

	#access geometry and topology information ( odb->rootAssembly->instances->(nodes, elements) )
	rootassembly = odb.rootAssembly
	instance = rootassembly.instances
	#access attribute information
	step = odb.steps
	#get instance & step information : Quantity and all names
	allinstancestr = str(instance)
	autoins = allinstancestr.split("'")
	inslen = len(autoins)/4
	instance_N = range(0,inslen)
	allstepstr = str(step)
	autostep = allstepstr.split("'")
	steplen = len(autostep)/4
	step_N = range(0,steplen)
	
	for i in input_step:
		if(steplen < int(i)):
			print "input step exceeds the range of steps"
			os._exit(0)
	for i in input_instance:
		if(inslen < int(i)):
			print "input instance exceeds the range of instances"
			os._exit(0)
		

	#step cycle
	for step_i in input_step:
		n = int(step_i)*4+1
		stepname = autostep[n]
		print "Step: ",stepname
		#instance cycle
		for ins_i in input_instance:
			n = int(ins_i)*4+1
			instancename = autoins[n]
			print "Instance: ",instancename
			
			#access nodes & elements
			node = instance[instancename].nodes
			element = instance[instancename].elements
			n_nodes = len(node)
			n_elements = len(element)
			#access attribute(fieldOutputs) information
			frame = step[stepname].frames
			
			
			#compute the number of element of each block
			p_elements = n_elements/piecenum + 1
			lp_elements = n_elements - (p_elements*(piecenum-1))  #last block
			
			
			#match nodes' label and its order in sequence (for empty nodes in tetra mesh)
			MLN = node[n_nodes-1].label
			TOTAL=[]
			#read node in sequence, and get the largest label of node(non-empty) 
			#MLN is the max label of nodeset
			for i in node:
				TOTAL.append(i.label)
				if(i.label > MLN):
					MLN = i.label
			#match (the key)
			L=[]
			n = 0
			for i in range(MLN): 
				L.append(0)
			for i in TOTAL:
				L[i-1] = n
				n += 1

			
			#frame cycle
			for i_frame in input_frame:
				
				#Detect whether the input frame is out of range
				try:
					TRY = odb.steps[stepname].frames[int(i_frame)]
				except:
					print "input frame exceeds the range of frames" 
					os._exit(0)
					break
				
				#Access a frame
				N_Frame = odb.steps[stepname].frames[int(i_frame)]
				print "Frame:",i_frame
				
				#create array for store result data temporarily
				# Vector-U,A,V,RF 
				L0=[] 
				# Tensors-S
				L1=[]
				# Tensors-LE
				L2=[]
				# Tensors-PE
				L3=[]
				# Scalars-PEEQ
				L4=[]
				for i in range(MLN): 
					L0.append([0,0,0,0,0,0,0,0,0,0,0,0])
					L1.append([0,0,0,0,0,0,0,0,0,0,0,0,0,0])
					L2.append([0,0,0,0,0,0,0,0,0])
					L3.append([0,0,0,0,0,0,0,0,0])
					L4.append([0,0])
				
				
				print "Reading U, A, V, RF ......"
				time1 = time()
				#Access Spatial displacement
				displacement = N_Frame.fieldOutputs['U']
				fieldValues = displacement.values
				for valueX in fieldValues :
					i = valueX.nodeLabel
					L0[i-1][0] = valueX.data[0]
					L0[i-1][1] = valueX.data[1]
					L0[i-1][2] = valueX.data[2]
				#Access Spatial acceleration
				acceleration = N_Frame.fieldOutputs['A']
				fieldValues = acceleration.values
				for valueX in fieldValues :
					i = valueX.nodeLabel
					L0[i-1][3] = valueX.data[0]
					L0[i-1][4] = valueX.data[1]
					L0[i-1][5] = valueX.data[2]
				#Access Spatial velocity
				velocity = N_Frame.fieldOutputs['V']
				fieldValues = velocity.values
				for valueX in fieldValues :
					i = valueX.nodeLabel
					L0[i-1][6] = valueX.data[0]
					L0[i-1][7] = valueX.data[1]
					L0[i-1][8] = valueX.data[2]
				#Access Reaction force
				Reaction_force = N_Frame.fieldOutputs['RF']
				fieldValues = Reaction_force.values	
				for valueX in fieldValues :
					i = valueX.nodeLabel
					L0[i-1][9] = valueX.data[0]
					L0[i-1][10] = valueX.data[1]
					L0[i-1][11] = valueX.data[2]	
				print "Time elapsed: ", time() - time1, "s"
				print "Reading Stress ......"
				time1 = time()
				#access Stress components
				Stress = N_Frame.fieldOutputs['S']
				node_Stress = Stress.getSubset(position=ELEMENT_NODAL)
				fieldValues = node_Stress.values
				for valueX in fieldValues :
					L1[valueX.nodeLabel-1][0] += 1
					L1[valueX.nodeLabel-1][1] += valueX.data[0]
					L1[valueX.nodeLabel-1][2] += valueX.data[1]
					L1[valueX.nodeLabel-1][3] += valueX.data[2]
					L1[valueX.nodeLabel-1][4] += valueX.data[3]
					L1[valueX.nodeLabel-1][5] += valueX.data[4]
					L1[valueX.nodeLabel-1][6] += valueX.data[5]
					L1[valueX.nodeLabel-1][7] += valueX.mises
					L1[valueX.nodeLabel-1][8] += valueX.maxPrincipal
					L1[valueX.nodeLabel-1][9] += valueX.midPrincipal
					L1[valueX.nodeLabel-1][10] += valueX.minPrincipal
					L1[valueX.nodeLabel-1][11] += valueX.press
					L1[valueX.nodeLabel-1][12] += valueX.tresca
					L1[valueX.nodeLabel-1][13] += valueX.inv3
				# can first ave
				print "Time elapsed: ", time() - time1, "s"
				print "Reading Logarithmic strain ......"
				time1 = time()
				#Logarithmic strain components
				Logarithmic_strain = N_Frame.fieldOutputs['LE']
				node_Logarithmic_strain = Logarithmic_strain.getSubset(position=ELEMENT_NODAL)
				fieldValues = node_Logarithmic_strain.values
				for valueX in fieldValues :
					L2[valueX.nodeLabel-1][0] += 1
					L2[valueX.nodeLabel-1][1] += valueX.data[0]
					L2[valueX.nodeLabel-1][2] += valueX.data[1]
					L2[valueX.nodeLabel-1][3] += valueX.data[2]
					L2[valueX.nodeLabel-1][4] += valueX.data[3]
					L2[valueX.nodeLabel-1][5] += valueX.data[4]
					L2[valueX.nodeLabel-1][6] += valueX.data[5]
					L2[valueX.nodeLabel-1][7] += valueX.maxPrincipal
					L2[valueX.nodeLabel-1][8] += valueX.minPrincipal
				print "Time elapsed: ", time() - time1, "s"
				print "Reading Plastic strain ......"
				time1 = time()
				#Plastic strain components
				Plastic_strain = N_Frame.fieldOutputs['PE']
				node_Plastic_strain = Plastic_strain.getSubset(position=ELEMENT_NODAL)
				fieldValues = node_Plastic_strain.values	
				for valueX in fieldValues :
					L3[valueX.nodeLabel-1][0] += 1
					L3[valueX.nodeLabel-1][1] += valueX.data[0]
					L3[valueX.nodeLabel-1][2] += valueX.data[1]
					L3[valueX.nodeLabel-1][3] += valueX.data[2]
					L3[valueX.nodeLabel-1][4] += valueX.data[3]
					L3[valueX.nodeLabel-1][5] += valueX.data[4]
					L3[valueX.nodeLabel-1][6] += valueX.data[5]
					L3[valueX.nodeLabel-1][7] += valueX.maxPrincipal
					L3[valueX.nodeLabel-1][8] += valueX.minPrincipal
				print "Time elapsed: ", time() - time1, "s"
				print "Reading Equivalent plastic strain ......"
				time1 = time()
				#Equivalent plastic strain
				Equivalent_plastic_strain = N_Frame.fieldOutputs['PEEQ']
				node_Equivalent_plastic_strain = Equivalent_plastic_strain.getSubset(position=ELEMENT_NODAL)
				fieldValues = node_Equivalent_plastic_strain.values
				for valueX in fieldValues :
					L4[valueX.nodeLabel-1][0] += 1
					L4[valueX.nodeLabel-1][1] += valueX.data
				print "Time elapsed: ", time() - time1, "s"	
				
				
				'''============================================================'''
				
				print "Partitionning model and writing vtk files ......"
				#piece cycle, to partion the model and create each piece for vtk files		
				for pn in range(piecenum):
					time1 = time()
					print "frame:",i_frame,"; block:",pn
					#Reorganization
					#Control&Storage
					#estimate whether the node has already existed
					stg_p = []
					#store the reorganized node for element
					stg_e = []
					#store the reorganized node for node
					stg_n = []
					for i in range(MLN):
						stg_p.append(-1)
					nodecount = 0
					#reorganize the node and element (reconstruct the mesh)
					if(pn == piecenum-1):
						M = range(pn*p_elements,n_elements)
					else:
						M = range(pn*p_elements,(pn+1)*p_elements)
					for i in M:
						for j in range(mesh_conner):
							k = element[i].connectivity[j] - 1
							if(stg_p[k] < 0): 
								stg_p[k] = nodecount
								stg_n.append(L[k]) 
								stg_e.append(nodecount)
								nodecount += 1
							else:
								stg_e.append(stg_p[k])
					#compute point quantity
					n_reop = len(stg_n)
					reop_N = range(0,len(stg_n))


					#create and open a VTK(.vtu) files
					if(piecenum > 1):
						outfile = open (os.path.join(vtk_path,odbname)+'_'+stepname+'_'+instancename+'f%03d'%int(i_frame)+' '+'p'+str(pn)+'.vtu','w')
					if(piecenum == 1):
						outfile = open (os.path.join(vtk_path,odbname)+'_'+stepname+'_'+instancename+'f%03d'%int(i_frame)+'.vtu','w')
					
					#<VTKFile>, including the type of mesh, version, and byte_order
					outfile.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'+'\n')
					#<UnstructuredGrid>
					outfile.write('<UnstructuredGrid>'+'\n')
					#<Piece>, including the number of points and cells
					if(pn == piecenum-1):
						outfile.write('<Piece NumberOfPoints="'+str(n_reop)+'"'+' '+'NumberOfCells="'+str(lp_elements)+'">'+'\n')
					else:
						outfile.write('<Piece NumberOfPoints="'+str(n_reop)+'"'+' '+'NumberOfCells="'+str(p_elements)+'">'+'\n')

					
					print "Writing Nodes ......"
					#<Points> Write nodes into vtk files
					displacement = N_Frame.fieldOutputs['U']
					fieldValues = displacement.values
					outfile.write('<Points>'+'\n')
					outfile.write('<DataArray type="Float64" NumberOfComponents="3" format="ascii">'+'\n')
					for i in reop_N:
						nt = stg_n[i]
						k = node[stg_n[i]].label-1
						X,Y,Z = node[nt].coordinates[0]+L0[k][0],node[nt].coordinates[1]+L0[k][1],node[nt].coordinates[2]+L0[k][2]
						outfile.write(' '+'%11.8e'%X+'  '+'%11.8e'%Y+'  '+'%11.8e'%Z+'\n')			
					outfile.write('</DataArray>'+'\n')
					outfile.write('</Points>'+'\n')
					#</Points>


					print "Writing Results data ......"
					#<PointData> Write results data into vtk files
					outfile.write("<"+"PointData"+" "+"Tensors="+'"'+"Stress_Components,Plastic_strain_components"+'"'\
					+" "+"Vevtors="+'"'+"Spatial_displacement,Reaction_force"+'"'\
					+" "+"Scalars="+'"'+"Equivalent_plastic_strain,Stress_Mises,Stress_Max_Principal,Stress_Mid_Principal,Stress_Min_Principal,Stress_Pressure,Stress_Tresca,Stress_Third_Invariant,Plastic_strain_Max_Principal,Plastic_strain_Min_Principal"+'"'+">"+'\n')
					
					#Stress components, <DataArray>
					outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Components"+'"'+" "+"NumberOfComponents="+'"'+"9"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
					for i in reop_N:
						k = node[stg_n[i]].label-1
						XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ = L1[k][1]/L1[k][0],L1[k][4]/L1[k][0],L1[k][6]/L1[k][0],L1[k][4]/L1[k][0],L1[k][2]/L1[k][0],L1[k][5]/L1[k][0],L1[k][6]/L1[k][0],L1[k][5]/L1[k][0],L1[k][3]/L1[k][0]
						outfile.write('%11.8e'%XX+' '+'%11.8e'%XY+' '+'%11.8e'%XZ+' '+'%11.8e'%YX+' '+'%11.8e'%YY+' '+'%11.8e'%YZ+' '+'%11.8e'%ZX+' '+'%11.8e'%ZY+' '+'%11.8e'%ZZ+'\n')
					outfile.write("</DataArray>"+'\n')
					#</DataArray>

					#Logarithmic strain components, <DataArray>
					outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Logarithmic_strain_components"+'"'+" "+"NumberOfComponents="+'"'+"9"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
					for i in reop_N:
						k = node[stg_n[i]].label-1
						XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ = L2[k][1]/L2[k][0],L2[k][4]/L2[k][0],L2[k][6]/L2[k][0],L2[k][4]/L2[k][0],L2[k][2]/L2[k][0],L2[k][5]/L2[k][0],L2[k][6]/L2[k][0],L2[k][5]/L2[k][0],L2[k][3]/L2[k][0]
						outfile.write('%11.8e'%XX+' '+'%11.8e'%XY+' '+'%11.8e'%XZ+' '+'%11.8e'%YX+' '+'%11.8e'%YY+' '+'%11.8e'%YZ+' '+'%11.8e'%ZX+' '+'%11.8e'%ZY+' '+'%11.8e'%ZZ+'\n')
					outfile.write("</DataArray>"+'\n')
					#</DataArray>
					
					#Plastic strain components, <DataArray>
					outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Plastic_strain_components"+'"'+" "+"NumberOfComponents="+'"'+"9"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
					for i in reop_N:
						k = node[stg_n[i]].label-1
						XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ = L3[k][1]/L3[k][0],L3[k][4]/L3[k][0],L3[k][6]/L3[k][0],L3[k][4]/L3[k][0],L3[k][2]/L3[k][0],L3[k][5]/L3[k][0],L3[k][6]/L3[k][0],L3[k][5]/L3[k][0],L3[k][3]/L3[k][0]
						outfile.write('%11.8e'%XX+' '+'%11.8e'%XY+' '+'%11.8e'%XZ+' '+'%11.8e'%YX+' '+'%11.8e'%YY+' '+'%11.8e'%YZ+' '+'%11.8e'%ZX+' '+'%11.8e'%ZY+' '+'%11.8e'%ZZ+'\n')
					outfile.write("</DataArray>"+'\n')
					#</DataArray>

					#Spatial displacement, <DataArray>
					outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Spatial_displacement"+'"'+" "+"NumberOfComponents="+'"'+"3"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
					for i in reop_N:
						k = node[stg_n[i]].label-1
						X,Y,Z = L0[k][0],L0[k][1],L0[k][2]
						outfile.write('%11.8e'%X+' '+'%11.8e'%Y+' '+'%11.8e'%Z+'\n')
					outfile.write("</DataArray>"+'\n')
					#</DataArray>
		
					#Spatial acceleration, <DataArray>
					outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Spatial_acceleration"+'"'+" "+"NumberOfComponents="+'"'+"3"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
					for i in reop_N:
						k = node[stg_n[i]].label-1
						X,Y,Z = L0[k][3],L0[k][4],L0[k][5]
						outfile.write('%11.8e'%X+' '+'%11.8e'%Y+' '+'%11.8e'%Z+'\n')
					outfile.write("</DataArray>"+'\n')
					#</DataArray>

					#Spatial velocity, <DataArray>
					outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Spatial_velocity"+'"'+" "+"NumberOfComponents="+'"'+"3"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
					for i in reop_N:
						k = node[stg_n[i]].label-1
						X,Y,Z = L0[k][6],L0[k][7],L0[k][8]
						outfile.write('%11.8e'%X+' '+'%11.8e'%Y+' '+'%11.8e'%Z+'\n')
					outfile.write("</DataArray>"+'\n')	
					#</DataArray>
		
					#Reaction force
					outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Reaction_force"+'"'+" "+"NumberOfComponents="+'"'+"3"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
					for i in reop_N:
						k = node[stg_n[i]].label-1
						X,Y,Z = L0[k][9],L0[k][10],L0[k][11]
						outfile.write('%11.8e'%X+' '+'%11.8e'%Y+' '+'%11.8e'%Z+'\n')
					outfile.write("</DataArray>"+'\n')	
					#</DataArray>
					
					#Equivalent plastic strain, <DataArray>
					outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Equivalent_plastic_strain"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
					for i in reop_N:
						k = node[stg_n[i]].label-1
						X = L4[k][1]/L4[k][0]
						outfile.write('%11.8e'%X+'\n')
					outfile.write('</DataArray>'+'\n')
					#</DataArray>

					#Stress Mises, <DataArray>
					outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Mises"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
					for i in reop_N:
						k = node[stg_n[i]].label-1
						X = L1[k][7]/L1[k][0]
						outfile.write('%11.8e'%X+'\n')
					outfile.write('</DataArray>'+'\n')
					#</DataArray>

					#Stress Max.Principal, <DataArray>
					outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Max_Principal"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
					for i in reop_N:
						k = node[stg_n[i]].label-1
						X = L1[k][8]/L1[k][0]
						outfile.write('%11.8e'%X+'\n')
					outfile.write('</DataArray>'+'\n')	
					#</DataArray>
		
					#Stress Mid.Principal, <DataArray>
					outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Mid_Principal"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
					for i in reop_N:
						k = node[stg_n[i]].label-1
						X = L1[k][9]/L1[k][0]
						outfile.write('%11.8e'%X+'\n')
					outfile.write('</DataArray>'+'\n')	
					#</DataArray>
			
					#Stress Min.Principal, <DataArray>
					outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Min_Principal"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
					for i in reop_N:
						k = node[stg_n[i]].label-1
						X = L1[k][10]/L1[k][0]
						outfile.write('%11.8e'%X+'\n')
					outfile.write('</DataArray>'+'\n')	
					#</DataArray>
			
					#Stress Pressure, <DataArray>
					outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Pressure"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
					for i in reop_N:
						k = node[stg_n[i]].label-1
						X = L1[k][11]/L1[k][0]
						outfile.write('%11.8e'%X+'\n')
					outfile.write('</DataArray>'+'\n')	
					#</DataArray>

					#Stress Tresca, <DataArray>
					outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Tresca"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
					for i in reop_N:
						k = node[stg_n[i]].label-1
						X = L1[k][12]/L1[k][0]
						outfile.write('%11.8e'%X+'\n')
					outfile.write('</DataArray>'+'\n')	
					#</DataArray>
			
					#Stress Third_Invariant, <DataArray>
					outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Third_Invariant"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
					for i in reop_N:
						k = node[stg_n[i]].label-1
						X = L1[k][13]/L1[k][0]
						outfile.write('%11.8e'%X+'\n')
					outfile.write('</DataArray>'+'\n')	
					#</DataArray>
		
					#Logarithmic_strain_Max_Principal, <DataArray>
					outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Logarithmic_strain_Max_Principal"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
					for i in reop_N:
						k = node[stg_n[i]].label-1
						X = L2[k][7]/L2[k][0]
						outfile.write('%11.8e'%X+'\n')
					outfile.write('</DataArray>'+'\n')	
					#</DataArray>
		
					#Logarithmic strain Min.Principal, <DataArray>
					outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Logarithmic_strain_Min_Principal"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
					for i in reop_N:
						k = node[stg_n[i]].label-1
						X = L2[k][8]/L2[k][0]
						outfile.write('%11.8e'%X+'\n')
					outfile.write('</DataArray>'+'\n')	
					#</DataArray>'''
		
					#Plastic strain Max.Principal, <DataArray>
					outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Plastic_strain_Max_Principal"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
					for i in reop_N:
						k = node[stg_n[i]].label-1
						X = L3[k][7]/L3[k][0]
						outfile.write('%11.8e'%X+'\n')
					outfile.write('</DataArray>'+'\n')	
					#</DataArray>
		
					#Plastic strain Min.Principal, <DataArray>
					outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Plastic_strain_Min_Principal"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
					for i in reop_N:
						k = node[stg_n[i]].label-1
						X = L3[k][8]/L3[k][0]
						outfile.write('%11.8e'%X+'\n')
					outfile.write('</DataArray>'+'\n')	
					#</DataArray>
					
					outfile.write("</PointData>"+'\n')
					#</PointData>
					
					
					print "Writing Cells ......"
					#<Cells> Write cells into vtk files
					outfile.write('<Cells>'+'\n')
					#Connectivity
					outfile.write('<DataArray type="Int32" Name="connectivity" format="ascii">'+'\n')
					if (mesh_type == 12):
						for i in range(len(stg_e)/8):
							outfile.write(str(stg_e[i*8])+' '+str(stg_e[i*8+1])+' '+str(stg_e[i*8+2])+' '+str(stg_e[i*8+3])+' '+str(stg_e[i*8+4])+' '+str(stg_e[i*8+5])+' '+str(stg_e[i*8+6])+' '+str(stg_e[i*8+7])+'\n')
					if (mesh_type == 10):
						for i in range(len(stg_e)/4):
							outfile.write(str(stg_e[i*4])+' '+str(stg_e[i*4+1])+' '+str(stg_e[i*4+2])+' '+str(stg_e[i*4+3])+'\n')
					outfile.write('</DataArray>'+'\n')
					#Offsets
					outfile.write('<DataArray type="Int32" Name="offsets" format="ascii">'+'\n')
					for i in range(len(stg_e)/mesh_conner):
						outfile.write(str(i*mesh_conner+mesh_conner)+'\n')
					outfile.write('</DataArray>'+'\n')
					#Type
					outfile.write('<DataArray type="UInt8" Name="types" format="ascii">'+'\n')
					for i in range(len(stg_e)/mesh_conner):
						outfile.write(str(mesh_type)+'\n')
					outfile.write('</DataArray>'+'\n')
					outfile.write('</Cells>'+'\n')
					#</Cells>
		
		
					#</Piece>
					outfile.write('</Piece>'+'\n')
					#</UnstructuredGrid>
					outfile.write('</UnstructuredGrid>'+'\n')
					#</VTKFile>
					outfile.write('</VTKFile>'+'\n')
				
					outfile.close()
					print "Time elapsed: ", time() - time1, "s"	
				
				'''====================================================================='''
				print "Creating .pvtu file for frame ", i_frame," ......"
				#create .pvtu files for parallel visualization
				if ( piecenum > 1 ):
					outfile = open (os.path.join(vtk_path,odbname)+'_'+stepname+'_'+'f%03d'%int(i_frame)+'.pvtu','w')
					
					#write the basic information for .pvtu files
					outfile.write('<?xml version="1.0"?>'+'\n')
					outfile.write('<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">'+'\n')
					outfile.write("<PUnstructuredGrid GhostLevel="+'"'+str(piecenum)+'"'+">"+'\n')
					#pointdata
					outfile.write("<"+"PPointData"+" "+"Tensors="+'"'+"Stress_Components,Logarithmic_strain_components,Plastic_strain_components"+'"'\
					+" "+"Vevtors="+'"'+"Spatial_displacement,Spatial_acceleration,Spatial_velocity,Reaction_force"+'"'\
					+" "+"Scalars="+'"'+"Equivalent_plastic_strain,Stress_Mises,Stress_Max_Principal,Stress_Mid_Principal,Stress_Min_Principal,Stress_Pressure,Stress_Tresca,Stress_Third_Invariant,Logarithmic_strain_Max_Principal,Logarithmic_strain_Min_Principal,Plastic_strain_Max_Principal,Plastic_strain_Min_Principal"+'"'+">"+'\n')
					
					outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Components"+'"'+" "+"NumberOfComponents="+'"'+"9"+'"'+" "+"/>"+'\n')
					outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Logarithmic_strain_components"+'"'+" "+"NumberOfComponents="+'"'+"9"+'"'+" "+"/>"+'\n')
					outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Plastic_strain_components"+'"'+" "+"NumberOfComponents="+'"'+"9"+'"'+" "+"/>"+'\n')
					outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Spatial_displacement"+'"'+" "+"NumberOfComponents="+'"'+"3"+'"'+" "+"/>"+'\n')
					outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Spatial_acceleration"+'"'+" "+"NumberOfComponents="+'"'+"3"+'"'+" "+"/>"+'\n')
					outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Spatial_velocity"+'"'+" "+"NumberOfComponents="+'"'+"3"+'"'+" "+"/>"+'\n')
					outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Reaction_force"+'"'+" "+"NumberOfComponents="+'"'+"3"+'"'+" "+"/>"+'\n')
					outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Equivalent_plastic_strain"+'"'+" "+"/>"+'\n')
					outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Mises"+'"'+" "+"/>"+'\n')
					outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Max_Principal"+'"'+" "+"/>"+'\n')
					outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Mid_Principal"+'"'+" "+"/>"+'\n')
					outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Min_Principal"+'"'+" "+"/>"+'\n')
					outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Pressure"+'"'+" "+"/>"+'\n')
					outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Tresca"+'"'+" "+"/>"+'\n')
					outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_Third_Invariant"+'"'+" "+"/>"+'\n')
					outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Logarithmic_strain_Max_Principal"+'"'+" "+"/>"+'\n')
					outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Logarithmic_strain_Min_Principal"+'"'+" "+"/>"+'\n')
					outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Plastic_strain_Max_Principal"+'"'+" "+"/>"+'\n')
					outfile.write("<"+"PDataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Plastic_strain_Min_Principal"+'"'+" "+"/>"+'\n')
					outfile.write("</PPointData>"+'\n')
					
					#points
					outfile.write("<PPoints>"+'\n')
					outfile.write("<PDataArray type="+'"'+"Float64"+'"'+" "+"NumberOfComponents="+'"'+"3"+'"'+"/>"+'\n')
					outfile.write("</PPoints>"+'\n')
					
					#write the path of each piece for reading it through the .pvtu file 
					for pn in range(piecenum):
						outfile.write("<Piece Source="+'"'+odbname+'_'+stepname+'_'+instancename+'f%03d'%int(i_frame)+' '+'p'+str(pn)+'.vtu'+'"'+"/>"+'\n')
					
					outfile.write("</PUnstructuredGrid>"+'\n')
					outfile.write("</VTKFile>")

					outfile.close()	
			
			odb.close()

	print "Total time elapsed: ", time() - starttime, "s"
