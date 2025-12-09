import os
import wget 

def download_meshes(data_folder): 

    download_paths = {"armadillo" : "https://raw.githubusercontent.com/alecjacobson/common-3d-test-models/refs/heads/master/data/armadillo.obj",
                      "koala"     : "https://raw.githubusercontent.com/odedstein/reach-for-the-arcs-code/refs/heads/main/data/koala.obj",
                      "bunny"     : "https://raw.githubusercontent.com/alecjacobson/common-3d-test-models/refs/heads/master/data/stanford-bunny.obj",
                     }
    mesh_files     = {"armadillo" : "armadillo.obj", 
                      "koala"     : "koala.obj",
                      "bunny"     : "bunny.obj"}
    mesh_paths = {m:os.path.join(data_folder,f) for m,f in mesh_files.items()}

    # Download meshes if they're not available
    for mesh_name in mesh_paths.keys():
        if not os.path.exists(mesh_paths[mesh_name]):
            url = download_paths[mesh_name]
            filename = wget.download(url,out=mesh_paths[mesh_name])
            print()
            print(f"Downloaded {mesh_name} to {filename}")
        else:
            print(f"Found .obj file for {mesh_name}")

    return mesh_paths
