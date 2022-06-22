module msh_frm_binvox
    
    import Gmsh: gmsh
    export mesh_from_binvox

    function mesh_from_binvox(filename, direct_from_binvox=true)
        gmsh.initialize()
        gmsh.open(filename)
        entities = gmsh.model.getEntities()
        dim, tag = entities[end]

        #print info
        type = gmsh.model.getType(dim, tag)
        name = gmsh.model.getEntityName(dim, tag)
        if length(name) > 0
            name = name * " "
        end
        println("Entity " * name, entities[end], " of type ", type)

        # Extract the voxels/ linear hexahedrons
        #(type 5 element in gmsh)[https://gmsh.info/doc/texinfo/gmsh.html#Node-ordering]
        function get_elems(elem_id, dim, tag)
            elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, tag)
            for i in range(start=1, stop=length(elemTypes))
                if elemTypes[i] == elem_id
                    elemTags = elemTags[i]
                    elemNodeTags = elemNodeTags[i]
                    return elemTags, elemNodeTags
                end
            end
            error("elem_id $elem_id not found")
            return nothing
        end

        #get node data
        nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes(dim, tag)

        #get cell data
        hexahedrons = 5
        cell_id, cell_vertice_ids = get_elems(hexahedrons, dim, tag)

        # binvox generated msh file has some errors. Hence, the following workarounds. If generated directly using gmsh or some other means, set direct_from_binvox to false.
        if direct_from_binvox

            # # correction for the incorrect first tag
            nodeTags[1] = 0
            node_id = Int.(nodeTags) .+ 1
            
            # # cell corrections
            cell_id[1] = 0
            cell_vertice_ids[1] = 0
            cell_id = Int.(cell_id) .+ 1
            cell_vertice_ids = Int.(cell_vertice_ids) .+ 1

        else
            node_id = Int.(nodeTags);
            cell_id = Int.(cell_id);
        end

        #reshape data
        nodeCoords = reshape(nodeCoords, (3,:))
        cell_vertice_ids = reshape(cell_vertice_ids, (8,:))

        #round coords upto 3 decimal places
        nodeCoords = round.(nodeCoords, digits=3)

        gmsh.clear()
        gmsh.finalize()

        return cell_id, cell_vertice_ids, node_id, nodeCoords
    end
end
