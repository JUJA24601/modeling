from numpy import array, sort
import h5py
import copy
import sys

# Set size of computational box (I don't know)
origin = (0, 0, 0)
L = (1, 1, 1)

# Read commandline and open original files
num_files = len(sys.argv)-1
h5files = []    # The opened h5-files
groups = []     # The variable names, e.g., ("rhoe" or "rhoi")
filenames = []  # Names of h5-files
datas = []      # Name of data for steady state. For transient the datafiles will be named 0, 1, 2, ...
for filename in sys.argv[1:]:
    assert "_" in filename or "-" in filename, "Assuming filename is of type rhoi00000_xxx.h5 or rhoi00000-xxx.h5. No underscore or dash in "+filename
    
    filenames.append(filename)
    h5files.append(h5py.File(filename, "r+"))
    found_underscore = filename.find("_") if filename.find("_") >=0 else 1e8
    found_dash = filename.find("-") if filename.find("-") >=0 else 1e8
    if found_underscore < found_dash:
        groups.append(filename.split("_")[0].rstrip("0"))
#        groups.append("field")
    else:
        groups.append(filename.split("-")[0].rstrip("0"))
#        groups.append("field")
    datas.append(list(h5files[-1].keys())[0])

# Create a group and and rename all data items using integer time steps
# Only done first time around and only for transient calculations
if len(h5files[0]) > 1:
    print("Warning, regrouping transient h5files")
    for var, f in zip(groups, h5files):
        if not var in f:
            f.create_group(var)
            i = 0
            for key in f:    
                if not str(key) == var:
                    f.move(key, "/"+var+"/"+str(i))
                    i += 1

steady = False if isinstance(h5files[0][list(h5files[0].keys())[0]], h5py.Group) else True

# Create the light xdmf file
xdmffile = """<?xml version="1.0" encoding="utf-8"?>
<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.1">
  <Domain>
    <Grid Name="Structured Grid" GridType="Collection" CollectionType="Temporal">
"""
timeattr = """      <Time TimeType="List"><DataItem Format="XML" Dimensions="{1}"> {0} </DataItem></Time>"""

attribute3D_transient = """
        <Attribute Name="{0}" Center="Node">
          <DataItem Format="HDF" NumberType="Float" Precision="{4}" Dimensions="{1} {1} {1}">
            {2}:/{0}/{3:04d}
          </DataItem>
        </Attribute>"""
        
attribute3D_steady = """
        <Attribute Name="{0}" Center="Node">
          <DataItem Format="HDF" NumberType="Float" Precision="{4}" Dimensions="{1} {1} {1}">
            {2}:/{3}
          </DataItem>
        </Attribute>"""
              
xf3d = copy.copy(xdmffile)
#if not steady
timesteps = list(h5files[0][groups[0]].keys())
timesteps = sort(array(timesteps).astype(int))
xf3d += timeattr.format(str(timesteps)[1:-1], len(timesteps))
dims = h5files[0][groups[0]+"/0000"].shape
    
#else:
#    timesteps = [0]
#    dims = h5files[0][datas[0]].shape
    
for tstep in timesteps:
    xf3d += """
      <Grid GridType="Uniform">
        <Geometry Type="ORIGIN_DXDYDZ">
          <DataItem DataType="UInt" Dimensions="3" Format="XML" Precision="4">{0} {1} {2}</DataItem>
          <DataItem DataType="Float" Dimensions="3" Format="XML" Precision="4">{3} {4} {5}</DataItem>
        </Geometry>""".format(origin[0], origin[1], origin[2], L[0], L[1], L[2])

    xf3d += """
        <Topology Dimensions="{0} {1} {2}" Type="3DCoRectMesh"/>""".format(*dims)
    prec = 8 # double precision single is 4
#    if not steady:
    for var, filename in zip(groups, filenames):
        xf3d += attribute3D_transient.format(var, dims[0], filename, tstep, prec)
#    else:
#        for var, filename, data in zip(groups, filenames, datas):
#            xf3d += attribute3D_steady.format(var, dims[0], filename, data, prec)
        
    xf3d += """  
      </Grid>
"""
xf3d += """    
    </Grid>
  </Domain>
</Xdmf>  
"""
#f.close()
xf = open("_".join(groups) + ".xdmf", "w")
xf.write(xf3d)
xf.close()
