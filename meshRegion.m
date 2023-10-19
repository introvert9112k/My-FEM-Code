
function [node,element] = meshRegion(pt1, pt2, pt3, pt4, numx, numy, elemType) 

switch elemType

    case 'Q4'           % here we generate the mesh of Q4 elements
        nnx=numx+1; %no of nodes in x direction, if there are n elements then there will be n+1 nodes
        nny=numy+1; %no of nodes in y direction
        node=square_node_array(pt1,pt2,pt3,pt4,nnx,nny);
        inc_u=1;
        inc_v=nnx;
        node_pattern=[ 1 2 nnx+2 nnx+1 ];
        [element]=make_elem(node_pattern,numx,numy,inc_u,inc_v);

    case 'Q8'           % here we generate a mesh of Q9 elements
        nnx=numx+1;
        nny=numy+1;
        node=square_node_array(pt1,pt2,pt3,pt4,nnx,nny);
        inc_u=1;
        inc_v=nnx;
        node_pattern=[ 1 2 nnx+2 nnx+1 ];
        element=make_elem(node_pattern,numx,numy,inc_u,inc_v);
        [element,node]=q4totq8(element,node,numx,numy);

    otherwise
        error('For now, only Q4 and Q9 are supported by the mesh generator');
end

end % END OF FUNCTION meshRegion 

