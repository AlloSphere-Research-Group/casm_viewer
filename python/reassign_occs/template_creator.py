#This is a library of functions to create a template for allolib to use in casmviewer
#For script instructions search for START OF SCRIPT%%%
import numpy as np
import json
import poscar
import os

"""
return (g , x, y) a*x + b*y = gcd(x,y)
"""
def egcd(a,b):
    s1= np.sign(a)
    s2= np.sign(b)
    p1=1
    p2=0
    tp1=0
    tp2=1
    q=0
    i1=abs(a)
    i2=abs(b)
    while (i2!=0):
        q= i1 // i2
        i1=i1 % i2
        i1,i2=i2,i1
        p1= p1 - q * tp1
        p1,tp1=tp1,p1
        p2= p2 - q * tp2
        p2,tp2=tp2,p2
    p1 *=s1
    p2 *= s2
    return (p1*a + p2*b,p1,p2)



"""This function returns a elementary hermite transformation E
    with determinant 1 suchs that E * [ a ; b ] = [g ; 0 ]
"""
def elementary_hermite_op(a , b, i, j):
    tmat=np.identity(3).astype(int)
    tgcf,p1,p2 = egcd(a, b)
    tmat[i,i]=p1
    tmat[i,j]=p2
    if tgcf==0:
        return np.identity(3).astype(int)
    tmat[j,i]= -b // tgcf
    tmat[j,j]= a // tgcf
    return tmat


"""This function takes in a transformation matrix as a numpy array
   and gives back the three matrices of the smith normal form as a tuple
   (U, S, V)
""" 
def smith_normal_form(transfmat):
    U = np.identity(3).astype(int)
    V = np.identity(3).astype(int)
    S = transfmat
    tmat = U;
    for j in range(3):
        for i in range(j+1,3):
            if (S[i,j]==0):
                continue
            tmat = elementary_hermite_op(S[j,j],S[i,j],j,i)
            S = np.dot(tmat,S)
            U = np.dot(U,np.linalg.inv(tmat))
        for i in range(j+2,3):
            if(S[j,i]==0):
                continue
            tmat = elementary_hermite_op(S[j,j+1],S[j,i],j+1,i)
            S = np.dot(S,tmat.transpose())
            V = np.dot(np.linalg.inv(tmat.transpose()),V)

    while (np.count_nonzero(S - np.diag(np.diagonal(S)))):
        b=0
        for b in range(2):
            if S[b,b+1]:
                break
        if (S[b,b] < 0):
            S[b,:] = -S[b,:]
            U[:,b] = -U[:,b]
        if (S[b,b]):
            q = S[b,b+1] // S[b,b]
            #if mod gives negative then negative number?
            if (S[b,b+1] % S[b,b] < 0):
                q-=1
            tmat = np.identity(3).astype(int)
            tmat[b+1,b]=-q
            S=np.dot(S,tmat.transpose())
            V=np.dot(np.linalg.inv(tmat.transpose()),V)
        else:
            tmat = np.identity(3).astype(int)
            tmat[b,b]=0
            tmat[b,b+1]=1
            tmat[b+1,b+1]=0
            tmat[b+1,b]=1
            S=np.dot(S,tmat.transpose())
            V=np.dot(np.linalg.inv(tmat.transpose()),V)

        if not S[b,b+1]:
            continue

        tmat = elementary_hermite_op(S[b,b],S[b,b+1],b,b+1)
        S=np.dot(S,tmat.transpose())
        V=np.dot(np.linalg.inv(tmat.transpose()),V)
        for j in range(2):
            tmat = elementary_hermite_op(S[j,j],S[j+1,j],j,j+1)
            S=np.dot(tmat,S)
            U=np.dot(U,np.linalg.inv(tmat))

            if (j+2 >=3):
                continue
            
            tmat = elementary_hermite_op(S[j,j+1],S[j,j+2],j+1,j+2)
            S=np.dot(S,tmat.transpose())
            V=np.dot(np.linalg.inv(tmat.transpose()),V)

    for j in range(3):
        if S[j,j] < 0:
            for i in range(3):
                S[j,i]=-S[j,i]
                U[i,j]=-U[i,j]

    for i in range(3):
        for j in range(i+1,3):
            if (S[i,i] > S[j,j]):
                S[[i,j]] = S[[j,i]]
                S[:,[i,j]] = S[:,[j,i]]
                U[:,[i,j]] = U[:,[j,i]]
                V[[i,j]] = V[[j,i]]
#                S[i],S[j]=S[j],S[i]
#                S[:,i],S[:,j]=S[:,j],S[:,i]
#                U[:,i],U[:,j] = U[:,j],U[:,i]
#                V[i],V[j] = V[j],V[i]

#for(i = 0; i < 3; i++) {
#      if(S(i, i) == 0) continue;
#      for(j = i + 1; j < 3; j++) {
#        if(S(j, j) % S(i, i) == 0) continue;
#        //Replace S(i,i), S(j,j) by their gcd and lcm respectively.
#        tmat = DerivedOut::Identity();
#        DerivedOut tmat2(tmat);
#        Scalar a(S(i, i)), b(S(j, j)), c, d, tgcf;
#        tgcf = extended_gcf(a, b, c, d);
#        tmat(i, i) = 1;
#        tmat(i, j) = d;
#        tmat(j, i) = -b / tgcf;
#        tmat(j, j) = (a * c) / tgcf;
#
#        tmat2(i, i) = c;
#        tmat2(i, j) = 1;
#        tmat2(j, i) = -(b * d) / tgcf;
#        tmat2(j, j) = a / tgcf;
#        S = tmat * S * tmat2.transpose();
#        U = U * inverse(tmat);
#        V = inverse(tmat2.transpose()) * V;
#      }
#    }
    for i in range(3):
        if (S[i,i]==0):
            continue
        for j in range(i+1,3):
            if (S[j,j] % S[i,i] == 0):
                continue
            tmat = np.identity(3).astype(int)
            tmat2 = tmat.copy()
            tgcf,c,d = egcd(S[i,i],S[j,j])
            tmat[i,i]=1
            tmat[i,j]=d
            tmat[j,i]=-S[j,j] // tgcf
            tmat[j,j]=(S[i,i] * c) // tgcf

            
            tmat2[i,i]=c
            tmat2[i,j]=1
            tmat2[j,i]=-(S[j,j]* d) // tgcf
            tmat2[j,j]=S[i,i] // tgcf

            S = np.dot(np.dot(tmat,S),tmat2.transpose())
            U = np.dot(U,np.linalg.inv(tmat))
            V = np.dot(np.linalg.inv(tmat2.transpose()),V)
            
#    return (U.astype(int), S.astype(int), V.astype(int))
    return (np.rint(U).astype(int), np.rint(S).astype(int), np.rint(V).astype(int))

"""This function takes in the three components of the smith
   normal form decomposition of the transformation matrix
   as numpy arrays and returns a list of numpy arrays
   that represents the unitcells in sorted order
"""
def get_ordered_unit_cells(U, S, V):
    unitcells=[]
    transfmat=np.dot(np.dot(U,S),V)
    nvol=round(np.linalg.det(transfmat))
    Smat=S.copy();
    Smat[0,0]=nvol / Smat[0,0]
    Smat[1,1]=nvol / Smat[1,1]
    Smat[2,2]=nvol / Smat[2,2]
    for i in range(int(round(nvol))):
        mnp=np.array([ (i % (S[0,0]*S[1,1])) % S[0,0],
                        (i % (S[0,0]*S[1,1])) // S[0,0] ,
                        i // (S[0,0]*S[1,1])])
        ijk=np.dot(U,mnp.transpose())
        realijk=ijk
        vec2=np.dot(np.dot(np.dot(np.linalg.inv(V),Smat),np.linalg.inv(U)),ijk)
        vec2[0]=((vec2[0]% nvol)+nvol)% nvol
        vec2[1]=((vec2[1]% nvol)+nvol)% nvol
        vec2[2]=((vec2[2]% nvol)+nvol)% nvol
        realijk= (np.dot(transfmat , vec2)) // nvol
        unitcells.append(realijk)
    return unitcells

"""This functions takes in the dictionary that represents the prim.json file
   and the list of unitcells that represent the tiling of the supercell of interest
   as a list of numpy arrays and returns a list of numpy arrays that represents
   the cartesian coordinates of all the basis sites in the supercell.
"""
def get_ordered_scel_cart_coords(primjsondict,unitcells):
    basis=primjsondict["basis"]
    lat_row_mat=np.array(primjsondict["lattice_vectors"])
    coord_mode=primjsondict["coordinate_mode"]
    cart_coords=[]
    #l mod scel_volume (length of unitcells) indicates which unitcell to shift to
    #l / scel_volume (length of unitcells) indicates which basis site in prim we are tiling
    for l in range(len(basis)*len(unitcells)):
           if coord_mode=="Fractional":
               to_add=np.dot(lat_row_mat.transpose(),np.array(basis[l // len(unitcells)]["coordinate"]).transpose())
           else:
               to_add=np.array(basis[l // len(unitcells)]["coordinate"]).transpose()
           cart_coords.append(np.dot(lat_row_mat.transpose(),unitcells[l % len(unitcells)].transpose()) + to_add)
    return cart_coords

"""This function takes in the filenames of the prim.json and the transformation matrix
    and reads in the data and feeds in to the functions that give a template lattice
    and a list of ordered template cartesian coordinates of the atoms in the transformed supercell.
    if the files are malformed this will give undefined behavior.
"""
def create_template_from_files(primjsonfilename,transfmatfilename):
    primjson=dict()
    with open(primjsonfilename) as f:
        primjson=json.load(f)
        with open(transfmatfilename) as mat_file:
            transfmat=np.loadtxt(mat_file)
            U,S,V=smith_normal_form(transfmat)
            unitcells=get_ordered_unit_cells(U,S,V)
            superlattice=np.dot(np.array(primjson["lattice_vectors"]).transpose(),np.dot(np.dot(U,S),V)).transpose()
            supercart_coords=get_ordered_scel_cart_coords(primjson,unitcells)
    return (superlattice,supercart_coords)
    
# def write_template(superlattice, supercart_coords, out_name = "template_POSCAR"):
#     #create a dummy template poscar to read in and populate fields of a poscar
#     #using poscar.read
#     dummy_poscar_string="New structure\n"\
#     "1.0\n"\
#     "     1.0    0.0    0.0\n"\
#     "     0.0    1.0    0.0\n"\
#     "     0.0    0.0    1.0\n"\
#     "     B\n"\
#     "     1\n"\
#     "Direct\n"\
#     "     0.000000000         0.000000000         0.000000000\n"

#     dummyposfile = open("dummypos",'w')
#     dummyposfile.write(dummy_poscar_string)
#     dummyposfile.close()
#     p=poscar.Poscar("dummypos")
#     #overwrite the contents of dummypos
#     p._lattice=superlattice
#     p.coord_mode = 'Cartesian'
#     new_basis = [ poscar.Site(True,x,"","X") for x in supercart_coords ]
#     p.basis = new_basis
#     p.type_atoms = ["X"]
#     p.type_atoms_alias = ["X"]
#     p.num_atoms= [len(p.basis)]
#     #write the template to template_POSCAR
#     p.write(out_name)
    
#     os.remove('dummypos') # remove temporary file

def write_netcdf(superlattice, supercart_coords, out_name = 'template.nc'):
    try:
        import netCDF4
    except:
        print("netCDF4 for python not available. Not saving netcdf template")
        return

    ncfile = netCDF4.Dataset(out_name, mode='w', format='NETCDF4') 
    # create complex128 numpy structured data type
    atom = np.dtype([('x',np.float32),('y',np.float32),('z',np.float32)])
    atom_t = ncfile.createCompoundType(atom,'atom')
    index = ncfile.createDimension('index')
    atoms_var = ncfile.createVariable('atoms_var',atom_t,('index'))

    lattice_dim = ncfile.createDimension('superlattice', size=9)
    lattice_var = ncfile.createVariable('lattice_var', np.float32,('superlattice',))
    print(superlattice)
    lattice_var[:] = superlattice.flat

    atoms_var[len(supercart_coords) - 1] = (0,0,0) # Does this serve to force preallocate?
    for i, coord in enumerate(supercart_coords):
        atoms_var[i] = (coord[0], coord[1], coord[2])

#The script part begins here START OF SCRIPT%%%
#This script expects to be run in a directory that contains a file named prim.json and transfmat
#The form of prim.json is {
#  "basis" : [
#    {
#      "coordinate" : [ 0.417374007600, 0.417374007600, 0.747877977200 ],
#      "occupant_dof" : [ "Va", "Na" ]
#    },
#    {
#      "coordinate" : [ 0.750707340934, 0.750707340933, 0.747877977200 ],
#      "occupant_dof" : [ "Va", "Na" ]
#    },
#    {
#      "coordinate" : [ 0.584040674266, 0.584040674267, 0.247877977200 ],
#      "occupant_dof" : [ "Co" ]
#    },
#    {
#      "coordinate" : [ 0.979766310928, 0.979766310927, 0.060701067218 ],
#      "occupant_dof" : [ "O" ]
#    },
#    {
#      "coordinate" : [ 0.188315037606, 0.188315037606, 0.435054887182 ],
#      "occupant_dof" : [ "O" ]
#    }
#  ],
#  "coordinate_mode" : "Fractional",
#  "description" : "P3_NaxCoO2",
#  "lattice_vectors" : [
#    [ 2.888862656049, 0.000000000000, 0.000000000000 ],
#    [ 1.444431328024, 2.501828448183, 0.000000000000 ],
#    [ 1.444431328024, 0.833942816061, 5.281161912909 ]
#  ],
#  "title" : "P3_NaxCoO2"
#}
#The form of transfmat is 
#-1 3 2
#0 5 1
#-3 4 6
#matrix must be integers and have integer determinant

if __name__ == "__main__":
    #TincArgumentParser to allow command line args or json config as only command line flag
    import sys, os
    from tinc import *
    parser = TincArgumentParser(description='POSCAR Template generator')

    parser.add_argument('__input_dir', type=str, default="./", nargs='?')
    parser.add_argument('__input_names', type=str, default="cached_output/transfmat", nargs='?')
    parser.add_argument('__output_names', type=str, default="template_POSCAR", nargs='?')
    parser.add_argument('__output_dir', type=str, default="", nargs='?')

    args = parser.get_args()

    # Try different prim names if not provided
    
    prim_name = args['__input_dir'] + "prim.json"
    if not os.path.exists(prim_name):
        prim_name = args['__input_dir'] + "../prim.json"
    if not os.path.exists(prim_name):
        prim_name = args['__input_dir'] + "prim_labels.json"
    if not os.path.exists(prim_name):
        prim_name = args['__input_dir'] + "../prim_labels.json"
        
    print("Using prim: " + prim_name, file = sys.stderr)

    lat, coords = create_template_from_files(prim_name, args['__input_dir']  + args['__input_names'][0])
    
    if args["__verbose"]:
        print("Writing output:" + args['__output_dir'] + args['__output_names'][0])
#    write_template(lat,coords, args['__output_dir'] + args['__output_name'])

    write_netcdf(lat,coords, args['__output_dir'] +  args['__output_names'][0])

