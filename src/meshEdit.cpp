#include <float.h>
#include <assert.h>
#include "meshEdit.h"
#include "mutablePriorityQueue.h"

namespace CMU462 {

  VertexIter HalfedgeMesh::splitEdge(EdgeIter e0) {

    // TODO: (meshEdit)
    // This method should split the given edge and return an iterator to the
    // newly inserted vertex. The halfedge of this vertex should point along
    // the edge that was split, rather than the new edges.
	if(!e0->isBoundary())
   {
   		HalfedgeIter h = e0->halfedge();
		int count = 0;
		std::vector<HalfedgeIter> boundaries;
	
		do
		{
			h = h->next();
			count++;
			boundaries.push_back(h);
		}while(h->next() != e0->halfedge());
	
		h = e0->halfedge()->twin();
	
		do
		{
			h = h->next();
			count++;
			boundaries.push_back(h);
		}while(h->next() != e0->halfedge()->twin());
	
		if(count == 4)
		{
			VertexIter v = newVertex();
			v->position = e0->centroid();
		
			EdgeIter e1 = e0;
			EdgeIter e2 = newEdge();
			EdgeIter e3 = newEdge();
			EdgeIter e4 = newEdge();
		
			FaceIter f1 = e0->halfedge()->face();
			FaceIter f4 = e0->halfedge()->twin()->face(); 
			FaceIter f3 = newFace();
			FaceIter f2 = newFace();
			
			HalfedgeIter h1 = e0->halfedge();
			HalfedgeIter h1_t = e0->halfedge()->twin();
			HalfedgeIter h2 = newHalfedge();
			HalfedgeIter h2_t = newHalfedge();
			HalfedgeIter h3 = newHalfedge();
			HalfedgeIter h3_t = newHalfedge();
			HalfedgeIter h4 = newHalfedge();
			HalfedgeIter h4_t = newHalfedge();
			
			
			h1->next() = boundaries[0];
			h1->twin() = h1_t;
			h1->edge() = e1;
			h1->face() = f1;
			h1->vertex() = v;
			
			f1->halfedge() = h1;
			v->halfedge() = h1;
			
			h1_t->next() = h4;
			h1_t->twin() = h1;
			h1_t->face() = f4;
			h1_t->edge() = e1;
			h1_t->vertex() = boundaries[3]->twin()->vertex();
			
			f4->halfedge() = h1_t;;
			
			boundaries[0]->next() = h2_t;
			boundaries[0]->face() = f1;
			
			h2_t->next() = h1;
			h2_t->twin() = h2;
			h2_t->edge() = e2;
			h2_t->face() = f1;
			h2_t->vertex() = boundaries[0]->twin()->vertex();
			
			boundaries[0]->twin()->vertex()->halfedge() = h2_t;
			
			h2->next() = boundaries[1];
			h2->twin() = h2_t;
			h2->edge() = e2;
			h2->face() = f2;
			h2->vertex() = v;
			
			e2->halfedge() = h2;
			f2->halfedge() = h2;
			
			boundaries[1]->next() = h3_t;
			boundaries[1]->face() = f2;
			
			h3_t->next() = h2;
			h3_t->twin() = h3;
			h3_t->edge() = e3;
			h3_t->face() = f2;
			h3_t->vertex() = boundaries[1]->twin()->vertex();
			
			boundaries[1]->twin()->vertex()->halfedge() = h3_t;
			
			h3->next() = boundaries[2];
			h3->twin() = h3_t;
			h3->edge() = e3;
			h3->face() = f3;
			h3->vertex() = v;
			
			e3->halfedge() = h3;
			f3->halfedge() = h3;
			
			boundaries[2]->next() = h4_t;
			boundaries[2]->face() = f3;
			
			h4_t->next() = h3;
			h4_t->twin() = h4;
			h4_t->edge() = e4;
			h4_t->face() = f3;
			h4_t->vertex() = boundaries[2]->twin()->vertex();
			
			boundaries[2]->twin()->vertex()->halfedge() = h4_t;
			
			h4->next() = boundaries[3];
			h4->twin() = h4_t;
			h4->edge() = e4;
			h4->face() = f4;
			h4->vertex() = v;
			
			e4->halfedge() = h4;
			f4->halfedge() = h4;
			
			boundaries[3]->next() = h1_t;
			boundaries[3]->face() = f4;
			
			boundaries[3]->twin()->vertex()->halfedge() = h1_t;
			
			return v;
		
    	
    	
		}else
		{
			return e0->halfedge()->vertex();	
		}
	
   }else
   {
   		HalfedgeIter h = e0->halfedge();
   		
   		if(h->isBoundary())
   		{
   			h = e0->halfedge()->twin();
   		}
		int count = 0;
		std::vector<HalfedgeIter> boundaries;
		
		HalfedgeIter prev_he = h->twin();
		HalfedgeIter next_he = h->twin()->next();
		
		do{
			prev_he = prev_he->next();
		}while(prev_he->next() != h->twin());
	
		HalfedgeIter h_sudo = h;
		do
		{
			h = h->next();
			count++;
			boundaries.push_back(h);
		}while(h->next() != h_sudo);
		
		h = h->next();
		
		if(count == 2)
		{
			VertexIter v = newVertex();
			v->position = e0->centroid();
			
			EdgeIter e1 = e0;
			EdgeIter e2 = newEdge();
			EdgeIter e3 = newEdge();
			
			FaceIter f1 = e0->halfedge()->face();
			FaceIter f2 = newFace();
			
			HalfedgeIter h1 = e0->halfedge();
			HalfedgeIter h1_t = e0->halfedge()->twin();
			
			HalfedgeIter h2 = newHalfedge();
			HalfedgeIter h2_t = newHalfedge();
			
			HalfedgeIter h3 = newHalfedge();
			HalfedgeIter h3_t = newHalfedge();
			
			v->halfedge() = h1;
			
			h1->next() = boundaries[0];
			h1->twin() = h1_t;
			h1->edge() = e1;
			h1->face() = f1;
			h1->vertex() = v;
			
			h1_t->next() = h3_t;
			h1_t->twin() = h1;
			h1_t->face() = prev_he->face();
			h1_t->edge() = e1;
			h1_t->vertex() = boundaries[0]->vertex();
			
			boundaries[0]->vertex()->halfedge() = h1_t;
			
			prev_he->next() = h1_t;
			
			e1->halfedge() = h1;
			f1->halfedge() = h1;
			
			h2->next() = boundaries[1];
			h2->twin() = h2_t;
			h2->edge() = e2;
			h2->face() = f2;
			h2->vertex() = v;
			
			h2_t->next() = h1;
			h2_t->twin() = h2;
			h2_t->edge() = e2;
			h2_t->face() = f1;
			h2_t->vertex() = boundaries[0]->twin()->vertex();
			
			e2->halfedge() = h2;
			f2->halfedge() = h2;
			
			boundaries[0]->next() = h2_t;
			boundaries[0]->face() = f1;
			
			h3->next() = h2;
			h3->twin() = h3_t;
			h3->edge() = e3;
			h3->face() = f2;
			h3->vertex() = boundaries[1]->twin()->vertex();
			
			boundaries[1]->twin()->vertex()->halfedge() = h3;
			
			h3_t->next() = next_he;
			h3_t->twin() = h3;
			h3_t->edge() = e3;
			h3_t->face() = next_he->face();
			h3_t->vertex() = v;
			
			boundaries[1]->next() = h3;
			boundaries[1]->face() = f2;
			
			
			return v;
			
		}else
		{
			return e0->halfedge()->vertex();		
		}
	
   	
   }
   
  }

  VertexIter HalfedgeMesh::collapseEdge(EdgeIter e) {

    // TODO: (meshEdit)
    // This method should collapse the given edge and return an iterator to
    // the new vertex created by the collapse.
    std::vector<HalfedgeIter> neighbours1;
    std::vector<HalfedgeIter> neighbours2;
    
    HalfedgeIter h1 = e->halfedge()->next();
    HalfedgeIter h2 = e->halfedge()->twin()->next();
        
    do
    {
    	neighbours1.push_back( h1 );
    	h1 = h1->twin()->next();
    }while(h1 != e->halfedge()->twin());
    
    do
    {
    	neighbours2.push_back( h2 );
    	h2 = h2->twin()->next();
    }while(h2 != e->halfedge());
    
    if(neighbours1.size() >= 2 || neighbours2.size() >= 2)
    {
    if( neighbours1.back()->twin()->vertex() == neighbours2.front()->twin()->vertex() && neighbours1.size() >= 2) 
    {

		cout<<"size1: "<<neighbours1.size()<<endl;
		HalfedgeIter h1_prev = neighbours1[neighbours1.size() - 2]->twin();
		HalfedgeIter h2_prev = neighbours2.front();
		
// 		printf("yellllllllllllllllllllllllloooooooooo\n");
		//halfedge assignment
		h1_prev->next() = e->halfedge()->twin();
		h2_prev->next() = neighbours1.back()->next();
		
		//face assignment
		e->halfedge()->twin()->face() = h1_prev->face();//neighbours1.back()->face();
		neighbours2.front()->face() = h1_prev->face();//neighbours1.back()->face();
		
		//neighbours1.back()->face()->halfedge() = neighbours1.back()->next();
		h1_prev->face()->halfedge() = h1_prev;
		
		//vertex assignment
		h1_prev->twin()->vertex()->halfedge() = h1_prev->twin();
		h2_prev->twin()->vertex()->halfedge() = h2_prev->twin();
// 		neighbours1.back()->vertex()->halfedge() = neighbours1[neighbours1.size() - 2];
// 		neighbours1.back()->twin()->vertex()->halfedge() =neighbours2.front()->twin();
		
		HalfedgeIter hTBD  = neighbours1.back();
		neighbours1.pop_back();
		
		deleteFace(hTBD->twin()->face());
		deleteEdge(hTBD->edge());
		deleteHalfedge(hTBD->twin());
		deleteHalfedge(hTBD);
		
		
    }
    
    if( neighbours1.front()->twin()->vertex() == neighbours2.back()->twin()->vertex() && neighbours2.size() >= 2)
    {
    	cout<<"size2: "<<neighbours2.size()<<endl;
		HalfedgeIter h2_prev = neighbours2[neighbours2.size() - 2]->twin();
		HalfedgeIter h1_prev = neighbours1.front();
		
		//halfedge assignment
		h2_prev->next() = e->halfedge();
		h1_prev->next() = neighbours2.back()->next();
		
		//face assignment
		e->halfedge()->face() = h2_prev->face();//neighbours2.back()->face();
		neighbours1.front()->face() = h2_prev->face();//neighbours2.back()->face();
		
		//neighbours2.back()->face()->halfedge() = neighbours2.back()->next();
		h2_prev->face()->halfedge() = h2_prev;
		
		//vertex assignment
		h1_prev->twin()->vertex()->halfedge() = h1_prev->twin();
		h2_prev->twin()->vertex()->halfedge() = h2_prev->twin();
// 		neighbours2.back()->vertex()->halfedge() = neighbours2[neighbours2.size() - 2];
// 		neighbours2.back()->twin()->vertex()->halfedge() =neighbours1.front()->twin();
		
		HalfedgeIter hTBD  = neighbours2.back();
		neighbours2.pop_back();
		
		deleteFace(hTBD->twin()->face());
		deleteEdge(hTBD->edge());
		deleteHalfedge(hTBD->twin());
		deleteHalfedge(hTBD);
		
		printf("here\n");
    }
    
    //edge deletion
	e->halfedge()->vertex()->position = e->centroid();
	
	VertexIter newVertex = e->halfedge()->vertex();
	VertexIter vertexTBD = e->halfedge()->twin()->vertex();
	
	neighbours2.back()->twin()->next() = neighbours1.front(); 
	neighbours1.back()->twin()->next() = neighbours2.front();

	newVertex->halfedge() = neighbours2.back();

	neighbours1.front()->face()->halfedge() = neighbours1.front();
	neighbours2.front()->face()->halfedge() = neighbours2.front();

	HalfedgeIter hTBD = e->halfedge();
	deleteEdge(hTBD->edge());
	deleteHalfedge(hTBD->twin());
	deleteHalfedge(hTBD);
	
	
	for(size_t i = 0; i < neighbours1.size(); i++)
	{
		neighbours1[i]->vertex() = newVertex;
	}
    
    deleteVertex(vertexTBD);
    return newVertex;
    
    }else
    {
    	return e->halfedge()->vertex();
    }
    
  }

  VertexIter HalfedgeMesh::collapseFace(FaceIter f) {

    // TODO: (meshEdit)
    // This method should collapse the given face and return an iterator to
    // the new vertex created by the collapse.

	if(!f->isBoundary() && f->degree() > 3)
    {
    	std::vector<HalfedgeIter> boundaries;
		HalfedgeIter h = f->halfedge();
	
		do{
			boundaries.push_back(h);
			h = h->next();
		}while(h != f->halfedge());

		VertexIter v = newVertex();
	
		v->position = f->centroid();
		
		v->halfedge() = boundaries[0]->twin()->next();
		
		
		for(size_t i = 0; i < boundaries.size(); i++)
		{
			boundaries[i]->twin()->next()->vertex() = v;
			
			i == boundaries.size() - 1 ? boundaries[0]->twin()->next()->twin()->next() = boundaries[i]->twin()->next() : boundaries[i + 1]->twin()->next()->twin()->next() = boundaries[i]->twin()->next(); 

			boundaries[i]->twin()->next()->face()->halfedge() = boundaries[i]->twin()->next();
		}
	
		for(size_t i = 0; i < boundaries.size(); i++)
		{
			deleteEdge( boundaries[i]->edge() );
			deleteVertex( boundaries[i]->vertex() );
			deleteHalfedge( boundaries[i]->twin() );
			deleteHalfedge( boundaries[i] );
		}
		
		deleteFace( f );
		return v;
		
    }else
    {
    	return f->halfedge()->vertex();
    }
    
  }

  FaceIter HalfedgeMesh::eraseVertex(VertexIter v) {

    // TODO: (meshEdit)
    // This method should replace the given vertex and all its neighboring
    // edges and faces with a single face, returning the new face.
	if(!v->isBoundary())
	{
		std::vector<VertexIter> neighbours;
		std::vector<HalfedgeIter> hedgeTBD;
		std::vector<HalfedgeIter> newhalfedge;
	
		HalfedgeIter h = v->halfedge();
	
		hedgeTBD.push_back(h->twin());
	
		VertexIter v1 = h->twin()->vertex();
	
		neighbours.push_back(v1);
	
		h = h->next();
	
		newhalfedge.push_back(h);
	
		while(!(h->twin()->vertex()->position == v1->position))//collecting all the neighbours
		{
			VertexIter v_ = h->twin()->vertex();
		
			if(v_->position == v->position)
			{
				newhalfedge.pop_back();
				hedgeTBD.push_back(h);
				h = h->twin()->next();
				newhalfedge.push_back(h);
			}else
			{
				neighbours.push_back(v_);
				h = h->next();
				newhalfedge.push_back(h);
			}
		}
	
		FaceIter f = newFace();
		f->halfedge() = newhalfedge[0];

		//establishing new connectivity
		for(size_t i = 0; i < newhalfedge.size(); i++)
		{
			if( i == newhalfedge.size() - 1)
			{
				newhalfedge[i]->next() = newhalfedge[0];
				newhalfedge[i]->face() = f;
				newhalfedge[i]->vertex() =  neighbours[i];
				neighbours[i]->halfedge() = newhalfedge[i];
			
			}else
			{
				newhalfedge[i]->next() = newhalfedge[i+1];
				newhalfedge[i]->face() = f;
				newhalfedge[i]->vertex() =  neighbours[i];
				neighbours[i]->halfedge() = newhalfedge[i];
			}
		
		} 
	
		//deleting old connectivity
		for(size_t i = 0; i < hedgeTBD.size(); i ++)
		{	
			deleteEdge( hedgeTBD[i]->edge() );
			deleteFace( hedgeTBD[i]->face() );
			deleteHalfedge( hedgeTBD[i]->twin() );
			deleteHalfedge( hedgeTBD[i] );
		}
	
		deleteVertex( v );
	
		return f;
	

	}else
	{
		if(v->halfedge()->face()->isBoundary())
		{
			v->halfedge()->twin()->face();
		}else
		{
			return v->halfedge()->face();
		}
	}
		
  }

  FaceIter HalfedgeMesh::eraseEdge( EdgeIter e ) {
    // TODO: (meshEdit)
    // This method should erase the given edge and return an iterator to the
    // merged face.
    
    if(e->isBoundary())
    {
    	printf("Its a boundary edge\n");
    	return e->halfedge()->face();
    }
    
    if(e->halfedge()->next() == e->halfedge() || e->halfedge()->twin() == e->halfedge() )
    {
    	printf("Cannot erase\n");
    	return e->halfedge()->face();
    }
    HalfedgeIter h1 = e->halfedge();
    HalfedgeIter h2 = h1->twin();
    
    if((h1->next() != h2 && h2->next() != h1) && h1->twin()->face() == h1->face())
	{
			printf("cannot delete this edge1\n");
			return e->halfedge()->face();
	}
		
    
    if((h1->next() == h2 || h2->next() == h1) && h1->face() == h1->twin()->face())
    {
    	HalfedgeIter h,h_prev;
    	if(h1->next() == h2)
    	{
    		h_prev = h1->next()->next();
    		h = h1->next()->next();
    		
    		do
			{
				h_prev = h_prev->next();

			}while(h_prev->next() != h1);
		
			h_prev->next() = h;
			h_prev->twin()->vertex()->halfedge() = h;
			h_prev->face()->halfedge() = h;
		
			deleteVertex(h1->twin()->vertex());
			deleteHalfedge(h1->twin());
			deleteHalfedge(h1);
			deleteEdge( e );
			return h->face();
			
		}else if(h2->next() == h1)
		{
			h_prev = h2->next()->next();
    		h = h2->next()->next();
    		
    		do
			{
				h_prev = h_prev->next();

			}while(h_prev->next() != h2);
		
			h_prev->next() = h;
			h_prev->twin()->vertex()->halfedge() = h;
			h_prev->face()->halfedge() = h;
		
			deleteVertex(h2->twin()->vertex());
			deleteHalfedge(h2->twin());
			deleteHalfedge(h2);
			deleteEdge( e );
			return h->face();
		}else
		{
			return e->halfedge()->face();
		}	
    	
    	
    }
    
	VertexIter v1 = h1->twin()->vertex();
	VertexIter v2 = h2->twin()->vertex();
	std::vector<HalfedgeIter> boundaries;
	std::vector<HalfedgeIter> twin_bound;
	
	
	do{
	
		h1 = h1->next();
		boundaries.push_back(h1);
		twin_bound.push_back(h1->twin());
	

	}while(!(h1->next()->twin()->vertex()->position == v1->position));
	
	
	do{
		h2 = h2->next();
		boundaries.push_back(h2);
		twin_bound.push_back(h2->twin());
		
	}while(!(h2->next()->twin()->vertex()->position == v2->position));
	
	
	FaceIter F = newFace();
	
	F->halfedge() = boundaries[0];
	
	
	for(size_t i = 0; i < boundaries.size(); i++)
	{
		if( i == boundaries.size() - 1)
		{
			
			boundaries[i]->next() = boundaries[0];

			boundaries[i]->face() =  F;
			boundaries[i]->vertex()->halfedge() = boundaries[i];
			boundaries[i]->edge()->halfedge() = boundaries[i];

		}else
		{

			boundaries[i]->next() = boundaries[i+1];

			boundaries[i]->edge()->halfedge() = boundaries[i];
			boundaries[i]->face() =  F;
			boundaries[i]->vertex()->halfedge() = boundaries[i]; 
		}
	}

	
	deleteFace( e->halfedge()->face() );
	deleteFace( e->halfedge()->twin()->face() );
	deleteHalfedge( e->halfedge()->twin() );
	deleteHalfedge( e->halfedge() );
	
	deleteEdge( e );
	return F;
	

  }

  EdgeIter HalfedgeMesh::flipEdge(EdgeIter e0) {

    // TODO: (meshEdit)
    // This method should flip the given edge and return an iterator to the
    // flipped edge.
    
     if(!e0->isBoundary())
    {
		HalfedgeIter h1 = e0->halfedge();
		HalfedgeIter h2 = e0->halfedge()->twin();
	
		std::vector<HalfedgeIter> boundary1;
		std::vector<HalfedgeIter> boundary2;
	
		do{
			h1 = h1->next();
			boundary1.push_back(h1);
		
		}while(h1->next() != e0->halfedge());
	
		do{
			h2 = h2->next();
			boundary2.push_back(h2);
		
		}while(h2->next() != e0->halfedge()->twin());//collecting boundary halfedges for both sides
	
		h1 = e0->halfedge();
		h2 = e0->halfedge()->twin();
		
		FaceIter f1 = h1->face();
		FaceIter f2 = h2->face();
		
		h1->next() = boundary1[1];
		h1->twin() = h2;
		h1->edge() = e0;
		h1->face() = f1;
		h1->vertex() = boundary2[0]->twin()->vertex();
		
		f1->halfedge() = h1;
		
		h2->next() = boundary2[1];
		h2->twin() = h1;
		h2->edge() = e0;
		h2->face() = f2;
		h2->vertex() = boundary1[0]->twin()->vertex();
		
		f2->halfedge() = h2;
		
		boundary1[0]->next() = h2;
		boundary1[0]->face() = f2;
		
		boundary1[boundary1.size() - 1]->next() = boundary2[0];
		boundary1[boundary1.size() - 1]->face() = f1;
		
		boundary2[0]->next() = h1;
		boundary2[0]->face() = f1;
		
		boundary2[boundary2.size() - 1]->next() = boundary1[0];
		boundary2[boundary2.size() - 1]->face() = f2;
		
		//e0->halfedge() = h1;
		
		return e0; 
	
    }else
    {
    	return e0;
    }

  }

  void HalfedgeMesh::subdivideQuad( bool useCatmullClark )
  {
    // Unlike the local mesh operations (like bevel or edge flip), we will perform
    // subdivision by splitting *all* faces into quads "simultaneously."  Rather
    // than operating directly on the halfedge data structure (which as you've seen
    // is quite difficult to maintain!) we are going to do something a bit nicer:
    //
    //    1. Create a raw list of vertex positions and faces (rather than a full-
    //       blown halfedge mesh).
    //
    //    2. Build a new halfedge mesh from these lists, replacing the old one.
    //
    // Sometimes rebuilding a data structure from scratch is simpler (and even more
    // efficient) than incrementally modifying the existing one.  These steps are
    // detailed below.

    // TODO Step I: Compute the vertex positions for the subdivided mesh.  Here we're
    // going to do something a little bit strange: since we will have one vertex in
    // the subdivided mesh for each vertex, edge, and face in the original mesh, we
    // can nicely store the new vertex *positions* as attributes on vertices, edges,
    // and faces of the original mesh.  These positions can then be conveniently
    // copied into the new, subdivided mesh.
    // [See subroutines for actual "TODO"s]
    if( useCatmullClark )
    {
      computeCatmullClarkPositions();
    }
    else
    {
      computeLinearSubdivisionPositions();
    }

    // TODO Step II: Assign a unique index (starting at 0) to each vertex, edge, and
    // face in the original mesh.  These indices will be the indices of the vertices
    // in the new (subdivided mesh).  They do not have to be assigned in any particular
    // order, so long as no index is shared by more than one mesh element, and the
    // total number of indices is equal to V+E+F, i.e., the total number of vertices
    // plus edges plus faces in the original mesh.  Basically we just need a one-to-one
    // mapping between original mesh elements and subdivided mesh vertices.
    // [See subroutine for actual "TODO"s]
    assignSubdivisionIndices();

    // TODO Step III: Build a list of quads in the new (subdivided) mesh, as tuples of
    // the element indices defined above.  In other words, each new quad should be of
    // the form (i,j,k,l), where i,j,k and l are four of the indices stored on our
    // original mesh elements.  Note that it is essential to get the orientation right
    // here: (i,j,k,l) is not the same as (l,k,j,i).  Indices of new faces should
    // circulate in the same direction as old faces (think about the right-hand rule).
    // [See subroutines for actual "TODO"s]
    vector< vector<Index> > subDFaces;
    vector< Vector3D > subDVertices;
    buildSubdivisionFaceList( subDFaces );
    buildSubdivisionVertexList( subDVertices );

    // TODO Step IV: Pass the list of vertices and quads to a routine that clears the
    // internal data for this halfedge mesh, and builds new halfedge data from scratch,
    // using the two lists.
    rebuild( subDFaces, subDVertices );
  }

  /**
   * Compute new vertex positions for a mesh that splits each polygon
   * into quads (by inserting a vertex at the face midpoint and each
   * of the edge midpoints).  The new vertex positions will be stored
   * in the members Vertex::newPosition, Edge::newPosition, and
   * Face::newPosition.  The values of the positions are based on
   * simple linear interpolation, e.g., the edge midpoints and face
   * centroids.
   */
  void HalfedgeMesh::computeLinearSubdivisionPositions()
  {
    // TODO For each vertex, assign Vertex::newPosition to
    // its original position, Vertex::position.

    // TODO For each edge, assign the midpoint of the two original
    // positions to Edge::newPosition.

    // TODO For each face, assign the centroid (i.e., arithmetic mean)
    // of the original vertex positions to Face::newPosition.  Note
    // that in general, NOT all faces will be triangles!
    
     for( VertexIter v = verticesBegin(); v != verticesEnd(); v++ )
    {
       v->newPosition = v->position;
    }
    
    
     for( EdgeIter e = edgesBegin(); e != edgesEnd(); e++ )
    {
       e->newPosition = e->centroid();
    }
    
	for( FaceIter f = facesBegin(); f != facesEnd(); f++ )
    {
       f->newPosition = f->centroid();
    }
    
    printf("calculated vertices\n");
    
  }

  /**
   * Compute new vertex positions for a mesh that splits each polygon
   * into quads (by inserting a vertex at the face midpoint and each
   * of the edge midpoints).  The new vertex positions will be stored
   * in the members Vertex::newPosition, Edge::newPosition, and
   * Face::newPosition.  The values of the positions are based on
   * the Catmull-Clark rules for subdivision.
   */
  void HalfedgeMesh::computeCatmullClarkPositions()
  {
    // TODO The implementation for this routine should be
    // a lot like HalfedgeMesh::computeLinearSubdivisionPositions(),
    // except that the calculation of the positions themsevles is
    // slightly more involved, using the Catmull-Clark subdivision
    // rules.  (These rules are outlined in the Developer Manual.)

    // TODO face
    
	for( FaceIter f = facesBegin(); f != facesEnd(); f++ )
    {
       f->newPosition = f->centroid();
    }
    
    for( EdgeIter e = edgesBegin(); e != edgesEnd(); e++ )
    {
       e->newPosition = (e->halfedge()->vertex()->position + e->halfedge()->twin()->vertex()->position + e->halfedge()->face()->newPosition + e->halfedge()->twin()->face()->newPosition)/4;
    }
    
    for( VertexIter v = verticesBegin(); v != verticesEnd(); v++ )
    {
    	Vector3D Q(0.f,0.f,0.f);
    	Vector3D R(0.f,0.f,0.f);
    	
    	HalfedgeIter h = v->halfedge();
    	
    	do{
    		Q = Q + h->face()->newPosition;
    		R = R + h->edge()->newPosition;
    		h = h->twin()->next();
    	}while(h != v->halfedge());	
    	
    	Q = Q / v->degree();
    	R = R / v->degree();
    	
    	v->newPosition = (Q + 2 * R + (v->degree() - 3) * v->position) / v->degree();
    }
    // TODO edges

    // TODO vertices
  }

  /**
   * Assign a unique integer index to each vertex, edge, and face in
   * the mesh, starting at 0 and incrementing by 1 for each element.
   * These indices will be used as the vertex indices for a mesh
   * subdivided using Catmull-Clark (or linear) subdivision.
   */
  void HalfedgeMesh::assignSubdivisionIndices()
  {
    // TODO Start a counter at zero; if you like, you can use the
    // "Index" type (defined in halfedgeMesh.h)

    // TODO Iterate over vertices, assigning values to Vertex::index

    // TODO Iterate over edges, assigning values to Edge::index

    // TODO Iterate over faces, assigning values to Face::index
    
    size_t count = 0;
    
	for( VertexIter v = verticesBegin(); v != verticesEnd(); v++ )
    {
       v->index = count;
       count++;
    }
    
    
    for( EdgeIter e = edgesBegin(); e != edgesEnd(); e++ )
    {
       e->index = count;
       count++;
    }
    
	for( FaceIter f = facesBegin(); f != facesEnd(); f++ )
    {
       f->index = count;
       count++;
    }
    
  }

  /**
   * Build a flat list containing all the vertex positions for a
   * Catmull-Clark (or linear) subdivison of this mesh.  The order of
   * vertex positions in this list must be identical to the order
   * of indices assigned to Vertex::newPosition, Edge::newPosition,
   * and Face::newPosition.
   */
  void HalfedgeMesh::buildSubdivisionVertexList( vector<Vector3D>& subDVertices )
  {
    // TODO Resize the vertex list so that it can hold all the vertices.
    
    // TODO Iterate over vertices, assigning Vertex::newPosition to the appropriate
    // location in the new vertex list.

    // TODO Iterate over edges, assigning Edge::newPosition to the appropriate
    // location in the new vertex list.

    // TODO Iterate over faces, assigning Face::newPosition to the appropriate
    // location in the new vertex list.
	
	for( VertexIter v = verticesBegin(); v != verticesEnd(); v++ )
    {
       subDVertices.push_back(v->newPosition);
    }
    
    
    for( EdgeIter e = edgesBegin(); e != edgesEnd(); e++ )
    {
       subDVertices.push_back(e->newPosition);
    }
    
	for( FaceIter f = facesBegin(); f != facesEnd(); f++ )
    {
       subDVertices.push_back(f->newPosition);
    }
    
    //printf("accumulated new vertices\n");
    
  }

  /**
   * Build a flat list containing all the quads in a Catmull-Clark
   * (or linear) subdivision of this mesh.  Each quad is specified
   * by a vector of four indices (i,j,k,l), which come from the
   * members Vertex::index, Edge::index, and Face::index.  Note that
   * the ordering of these indices is important because it determines
   * the orientation of the new quads; it is also important to avoid
   * "bowties."  For instance, (l,k,j,i) has the opposite orientation
   * of (i,j,k,l), and if (i,j,k,l) is a proper quad, then (i,k,j,l)
   * will look like a bowtie.
   */
  void HalfedgeMesh::buildSubdivisionFaceList( vector< vector<Index> >& subDFaces )
  {
    // TODO This routine is perhaps the most tricky step in the construction of
    // a subdivision mesh (second, perhaps, to computing the actual Catmull-Clark
    // vertex positions).  Basically what you want to do is iterate over faces,
    // then for each for each face, append N quads to the list (where N is the
    // degree of the face).  For this routine, it may be more convenient to simply
    // append quads to the end of the list (rather than allocating it ahead of
    // time), though YMMV.  You can of course iterate around a face by starting
    // with its first halfedge and following the "next" pointer until you get
    // back to the beginning.  The tricky part is making sure you grab the right
    // indices in the right order---remember that there are indices on vertices,
    // edges, AND faces of the original mesh.  All of these should get used.  Also
    // remember that you must have FOUR indices per face, since you are making a
    // QUAD mesh!

    // TODO iterate over faces
    // TODO loop around face
    // TODO build lists of four indices for each sub-quad
    // TODO append each list of four indices to face list
    
    for( FaceIter f = facesBegin(); f != facesEnd(); f++ )
    {
       if(!f->isBoundary())
       Size N = f->degree();
       HalfedgeIter h = f->halfedge();
       HalfedgeIter h_prev = h;
       
       do{
       		h_prev = h_prev->next();
       		
       }while(h_prev->next() != h);
       
       do{
       		std::vector<Index> quad( 4 );
       		quad[0] = h->vertex()->index;
       		quad[1] = h->edge()->index;
       		quad[2] = f->index;
       		quad[3] = h_prev->edge()->index;	
       		h = h->next();
       		h_prev = h_prev->next();
       		subDFaces.push_back(quad);
       }while(h != f->halfedge());
       
    }
    
    
  }

  void HalfedgeMesh::_bevel_fc_reposition_with_dist( vector<Vector3D>& orig, // list of vertex positions of the original face (before bevel)
                                                     vector<HalfedgeIter>& hs, // list of halfedges pointing from the vertices of the new, beveled face to the old, original face
                                                     double shift, // user-requested amount to shift the face in the normal direction
                                                     double inset ) // user-requested amount by which to inset (i.e., grow/shrink) the beveled face
  {
    // TODO Compute new vertex positions for the vertices of the beveled face.
    //
    // These vertices can be accessed via hs[i]->vertex()->position for i = 1, ..., hs.size()-1.
    //
    // The basic strategy here is to loop over the list of outgoing halfedges,
    // and use the preceding and next vertex position from the original mesh
    // (in the orig array) to compute an offset vertex position.
    //
    // Note that there is a 1-to-1 correspondence between halfedges in hs and vertex positions
    // in orig.  So, you can write loops of the form
    //
    // for( int i = 0; i < hs.size(); hs++ )
    // {
    //    Vector3D pi = orig[i]; // get the original vertex position correponding to vertex i
    // }
    //
    
    if(!hs[0]->twin()->next()->twin()->face()->isBoundary())
    {
    	size_t N = hs.size();
    
		Vector3D faceNormal = hs[0]->twin()->next()->twin()->face()->normal();
		faceNormal.unit();
		faceNormal = faceNormal*shift;
	
		for( size_t i = 0; i < hs.size(); i++)
		{
			size_t a = (i + N - 1) % N;
			size_t b = i;
			size_t c = (i+1) % N;
		
			Vector3D pa = orig[a];
			Vector3D pb = orig[b];
			Vector3D pc = orig[c];
		
			Vector3D newvect = (pa - pb) + (pc - pb);
			newvect.unit();
			newvect = newvect*inset;
		
			hs[i]->vertex()->position += newvect;
			hs[i]->vertex()->position -= faceNormal;
		
		} 
	
	}//
  }

  void HalfedgeMesh::_bevel_vtx_reposition_with_dist( Vector3D orig, // original vertex position, before the bevel
                                                      vector<HalfedgeIter>& hs, // list of halfedges pointing from the vertices of the new, beveled face to the neighbors of the original vertex
                                                      double inset ) // user-requested amount by which to inset (i.e., grow/shrink) the beveled face
  {
    // TODO Compute new vertex positions for the vertices of the beveled vertex.
    //
    // These vertices can be accessed via hs[i]->vertex()->position for i = 1, ..., hs.size()-1.
    //
    // The basic strategy here is to loop over the list of outgoing halfedges,
    // and use the preceding and next vertex position from the original mesh
    // (in the orig array) to compute an offset vertex position.
    //
    // Note that there is a 1-to-1 correspondence between halfedges in hs and vertex positions
    // in orig.  So, you can write loops of the form
    //
    // for( int i = 0; i < hs.size(); hs++ )
    // {
    //    Vector3D pi = orig[i]; // get the original vertex position correponding to vertex i
    // }
    //
    std::vector<Vector3D> vtx;
    
    for( size_t i = 0; i < hs.size(); i++ )
    {
    	vtx.push_back(hs[i]->twin()->vertex()->position);
    }
    
    size_t N = hs.size();
    
    for( size_t i = 0; i < hs.size(); i++)
    {
    	size_t a = (i + N - 1) % N;
    	size_t b = i;
    	size_t c = (i+1) % N;
    	
    	Vector3D newvect = hs[i]->vertex()->position - vtx[i];
    	newvect = newvect*inset;
    	
    	hs[i]->vertex()->position +=  newvect;
    	
    } 
  }

  void HalfedgeMesh::_bevel_edge_reposition_with_dist( vector<Vector3D>& origs,  // list of vertex positions of the neighbors of the two endpoints of the edge, before the bevel
                                                       vector<HalfedgeIter>& hs,  // list of halfedges pointing from the vertices of the new, beveled face to the neighbors of the endpoints of the old, original edge
                                                       double inset) // user-requested amount by which to inset (i.e., grow/shrink) the beveled face
  {
    // TODO Compute new vertex positions for the vertices of the beveled edge.
    //
    // These vertices can be accessed via hs[i]->vertex()->position for i = 1, ..., hs.size()-1.
    //
    // The basic strategy here is to loop over the list of outgoing halfedges,
    // and use the preceding and next vertex position from the original mesh
    // (in the orig array) to compute an offset vertex position.
    //
    // Note that there is a 1-to-1 correspondence between halfedges in hs and vertex positions
    // in orig.  So, you can write loops of the form
    //
    // for( int i = 0; i < hs.size(); hs++ )
    // {
    //    Vector3D pi = orig[i]; // get the original vertex position correponding to vertex i
    // }
    //
    
    size_t N = hs.size();
    
    for( size_t i = 0; i < hs.size(); i++)
    {
    	size_t a = (i + N - 1) % N;
    	size_t b = i;
    	size_t c = (i+1) % N;
    	
    	Vector3D newvect = hs[i]->vertex()->position - hs[i]->twin()->vertex()->position;
    	newvect = newvect*inset;
    	
    	hs[i]->vertex()->position -=  newvect;
    	
    }
  }

  FaceIter HalfedgeMesh::bevelVertex(VertexIter v) {

    // TODO This method should replace the vertex v with a face, corresponding to a bevel operation.
    // It should return the new face.  NOTE: This method is responsible for updating the *connectivity*
    // of the mesh only---it does not need to update the vertex positions.  These positions will be
    // updated in HalfedgeMesh::_bevel_vtx_reposition_with_dist (which you also have to implement!)
	
	HalfedgeIter h = v->halfedge();
	
	std::vector<HalfedgeIter> newhe;
	std::vector<HalfedgeIter> newhe_t;
	std::vector<HalfedgeIter> he;
	std::vector<HalfedgeIter> he_t;
	std::vector<EdgeIter> newedge;
	std::vector<VertexIter> newvtx;
	std::vector<FaceIter> face;
	
	FaceIter F = newFace();//new face
	
	do{
		
		he.push_back(h);
		he_t.push_back(h->twin());
		face.push_back(h->twin()->face());
		
		VertexIter v1 = newVertex();
		EdgeIter e1 = newEdge();
		HalfedgeIter h1 = newHalfedge();
		HalfedgeIter h1_t = newHalfedge();
		
		newvtx.push_back(v1);
		newedge.push_back(e1);
		newhe.push_back(h1);
		newhe_t.push_back(h1_t);
		
		h = h->twin()->next(); 
		
	}while(h != v->halfedge());
	
	for( size_t i = 0; i < he.size(); i++ ) 
	{
		he[i]->vertex() = newvtx[i];
		
	    he_t[i]->next() = newhe_t[i];
		
		i == he.size() - 1 ? newhe_t[i]->next() = he[0] : newhe_t[i]->next() = he[i + 1];
		newhe_t[i]->twin() = newhe[i];
		newhe_t[i]->edge() = newedge[i];
		newhe_t[i]->face() = face[i];
        newhe_t[i]->vertex() = newvtx[i];

		i == 0 ? newhe[i]->next() = newhe[he.size() - 1] : newhe[i]->next() = newhe[i - 1]; 		
		newhe[i]->twin() = newhe_t[i];
		newhe[i]->edge() = newedge[i];
		newhe[i]->face() = F;
		i == he.size() - 1 ? newhe[i]->vertex() = newvtx[0] : newhe[i]->vertex() = newvtx[i + 1];
		
		newedge[i]->halfedge() = newhe[i];
		newvtx[i]->halfedge() = he[i];
		newvtx[i]->position = v->position;
		face[i]->halfedge() = he_t[i];
		
	}
	
	F->halfedge() = newhe[0];
	deleteVertex(v);
	
    return F;
  }

  FaceIter HalfedgeMesh::bevelEdge(EdgeIter e) {

    // TODO This method should replace the edge e with a face, corresponding to a bevel operation.
    // It should return the new face.  NOTE: This method is responsible for updating the *connectivity*
    // of the mesh only---it does not need to update the vertex positions.  These positions will be
    // updated in HalfedgeMesh::_bevel_vtx_reposition_with_dist (which you also have to implement!)
	HalfedgeIter h1 = e->halfedge();
	HalfedgeIter h2 = e->halfedge()->twin();
	
	std::vector<HalfedgeIter> he;
	std::vector<HalfedgeIter> he_t;
	std::vector<FaceIter> face;
	std::vector<HalfedgeIter> newhe;
	std::vector<HalfedgeIter> newhe_t;
	std::vector<EdgeIter> newedge;
	std::vector<VertexIter> newvtx;
	std::vector<Vector3D> vtx;
	
	HalfedgeIter h = h1->next();
	do{
	
		he.push_back(h);
		he_t.push_back(h->twin());
		face.push_back(h->twin()->face());
		vtx.push_back(h->vertex()->position);
		
		VertexIter v = newVertex();
		EdgeIter e1 = newEdge();
		HalfedgeIter h_new = newHalfedge();
		HalfedgeIter h_new_t = newHalfedge();
		
		newhe.push_back(h_new);
		newhe_t.push_back(h_new_t);
		newedge.push_back(e1);
		newvtx.push_back(v);
	
		h = h->twin()->next();
		
	}while(h->edge() != e);
	
	h = h2;
	h = h->next();
	
	do{
	
	    he.push_back(h);
		he_t.push_back(h->twin());
		face.push_back(h->twin()->face());
		vtx.push_back(h->vertex()->position);
		
		VertexIter v = newVertex();
		EdgeIter e1 = newEdge();
		HalfedgeIter h_new = newHalfedge();
		HalfedgeIter h_new_t = newHalfedge();
		
		newhe.push_back(h_new);
		newhe_t.push_back(h_new_t);
		newedge.push_back(e1);
		newvtx.push_back(v);
	
		h = h->twin()->next();
		
	
	}while(h->edge() != e);
	
	FaceIter F = newFace();
	
	for( size_t i = 0 ; i < he.size(); i++ )
	{
		he[i]->vertex() = newvtx[i];
		
	    he_t[i]->next() = newhe_t[i];
		
		i == he.size() - 1 ? newhe_t[i]->next() = he[0] : newhe_t[i]->next() = he[i + 1];
		newhe_t[i]->twin() = newhe[i];
		newhe_t[i]->edge() = newedge[i];
		newhe_t[i]->face() = face[i];
        newhe_t[i]->vertex() = newvtx[i];

		i == 0 ? newhe[i]->next() = newhe[he.size() - 1] : newhe[i]->next() = newhe[i - 1]; 		
		newhe[i]->twin() = newhe_t[i];
		newhe[i]->edge() = newedge[i];
		newhe[i]->face() = F;
		i == he.size() - 1 ? newhe[i]->vertex() = newvtx[0] : newhe[i]->vertex() = newvtx[i + 1];
		
		newedge[i]->halfedge() = newhe[i];
		newvtx[i]->halfedge() = he[i];
		newvtx[i]->position = vtx[i];
		face[i]->halfedge() = he_t[i];
	}
	
	F->halfedge() = newhe[0];
	deleteVertex(e->halfedge()->vertex());
	deleteVertex(e->halfedge()->twin()->vertex());
	deleteHalfedge(e->halfedge()->twin());
	deleteHalfedge(e->halfedge());
	deleteEdge(e);
	
    return F;
  }

  FaceIter HalfedgeMesh::bevelFace(FaceIter f) {

    // TODO This method should replace the face f with an additional, inset face (and ring of faces around it),
    // corresponding to a bevel operation. It should return the new face.  NOTE: This method is responsible for
    // updating the *connectivity* of the mesh only---it does not need to update the vertex positions.  These
    // positions will be updated in HalfedgeMesh::_bevel_vtx_reposition_with_dist (which you also have to
    // implement!)
	
	if(!f->isBoundary())
	{
		HalfedgeIter h = f->halfedge();
		std::vector<HalfedgeIter> orig;
		std::vector<HalfedgeIter> newhefc;
		std::vector<HalfedgeIter> newhefc_t;
		std::vector<HalfedgeIter> newhe;
		std::vector<HalfedgeIter> newhe_t;
		std::vector<EdgeIter> newedg;
		std::vector<EdgeIter> newedg_d;
		std::vector<VertexIter> newvtx;
		std::vector<VertexIter> origvtx;
		std::vector<FaceIter> newfc;
	
		do
		{
			HalfedgeIter h1 = newHalfedge();
			HalfedgeIter h1_t = newHalfedge();
			HalfedgeIter h2 = newHalfedge();
			HalfedgeIter h2_t = newHalfedge();
			EdgeIter e1 = newEdge();
			EdgeIter e2 = newEdge();
			VertexIter v1 = newVertex();
			FaceIter f1 = newFace();
		
			newhefc.push_back(h1);
			newhefc_t.push_back(h1_t);
			newhe.push_back(h2);
			newhe_t.push_back(h2_t);
			newedg.push_back(e1);
			newedg_d.push_back(e2);
			newvtx.push_back(v1);
			newfc.push_back(f1);
			orig.push_back(h);
			origvtx.push_back(h->vertex());
		
			h = h->next();
		
		}while(h != f->halfedge());
	
		size_t N = orig.size();
	
		for(size_t i = 0; i < orig.size(); i++)
		{
			orig[i]->next() = newhe[i];
			orig[i]->face() = newfc[i];
		
			i == orig.size() - 1 ? newhefc[i]->next() = newhefc[0] : newhefc[i]->next() = newhefc[i+1];
			newhefc[i]->twin() = newhefc_t[i];
			newhefc[i]->face() = f;
			newhefc[i]->edge() = newedg[i];
			newhefc[i]->vertex() = newvtx[i];
		
			i == 0 ? newhefc_t[i]->next() = newhe_t[newhe.size() - 1] : newhefc_t[i]->next() = newhe_t[i - 1];
			newhefc_t[i]->twin() = newhefc[i];
			newhefc_t[i]->face() = newfc[i];
			newhefc_t[i]->edge() = newedg[i];
			i == orig.size() - 1 ?  newhefc_t[i]->vertex() = newvtx[0] : newhefc_t[i]->vertex() = newvtx[i+1];
			
			newhe[i]->next() = newhefc_t[i];
			newhe[i]->twin() = newhe_t[i];
			newhe[i]->face() = newfc[i];
			newhe[i]->edge() = newedg_d[i];
			newhe[i]->vertex() = orig[i]->twin()->vertex();
		
			i == orig.size() - 1 ? newhe_t[i]->next() = orig[0] : newhe_t[i]->next() = orig[i+1];
			newhe_t[i]->twin() = newhe[i];
			newhe_t[i]->edge() = newedg_d[i];
			i == orig.size() - 1 ? newhe_t[i]->face() = newfc[0] : newhe_t[i]->face() = newfc[i+1];
			i == orig.size() - 1 ? newhe_t[i]->vertex() = newvtx[0] : newhe_t[i]->vertex() = newvtx[i+1]; 
		
		
			size_t a = (i + N - 1) % N;
			size_t b = i;
			size_t c = (i+1) % N;
		
			Vector3D pa = origvtx[a]->position;
			Vector3D pb = origvtx[b]->position;
			Vector3D pc = origvtx[c]->position;
		
			Vector3D newvect = (pa - pb) + (pc - pb);
			newvect.unit();
			newvect = newvect * 0.01;
		
			newvtx[i]->halfedge() = newhefc[i];
			newvtx[i]->position = orig[i]->vertex()->position;//; + newvect;
			newedg[i]->halfedge() = newhefc[i];
			newedg_d[i]->halfedge() = newhe[i];
			newfc[i]->halfedge() = newhefc[i]->twin();	
		 
		}
	
		f->halfedge() = newhefc[0];
	
		return f;
	}
	
  }

  void HalfedgeMesh::splitPolygons(vector<FaceIter>& fcs) {
    for (auto f : fcs) splitPolygon(f);
  }

  void HalfedgeMesh::splitPolygon(FaceIter f) {
    // TODO triangulation
    
    if(!f->isBoundary())
    {
    	HalfedgeIter h = f->halfedge();
		std::vector<HalfedgeIter> helist;
	
		do{
			helist.push_back(h);
			h = h->next();
		}while(h != f->halfedge());//collecting all the halfedges
	
		if(helist.size() > 3)
		{
			int num_edges = helist.size() - 3;
			int list_size = helist.size();
			size_t list_cf = 0;
		
			while(num_edges > 0)
			{
				HalfedgeIter h1 = newHalfedge();
				HalfedgeIter h1_t = newHalfedge();
				EdgeIter e1 = newEdge();
				FaceIter f1 = newFace();
			
				HalfedgeIter h_prev  = helist[list_cf];
				
				do{
					h_prev = h_prev->next();
					//printf("d: %lu\n",helist.size());
	
				}while(h_prev->next()->next() != helist[list_cf]);
				
				//printf("I am here\n");
			
				HalfedgeIter h_next = h_prev->next();
				HalfedgeIter h_next_next = h_next->next();
			
			
				h1->next() = h_next;
				h1->twin() = h1_t;
				h1->face() = f1;
				h1->edge() = e1;
				h1->vertex() = h_next_next->twin()->vertex();//helist[list_cf]->vertex();
			
				f1->halfedge() = h1;
				e1->halfedge() = h1;
			
				h1_t->next() = h_next_next->next();
				h1_t->twin() = h1;
				h1_t->edge() = e1;
				h1_t->face() = f;
				h1_t->vertex() = h_next->vertex();
			
				f->halfedge() = h1_t;
			
				h_next_next->next() = h1;
				h_next_next->face() = f1;
			
				h_next->face() = f1;
				h_next->next() = h_next_next;
			
				h_prev->next() = h1_t;
				h_prev->face() = f;
			
				num_edges--;
				list_cf++;
	
				
						
				
			}
		}
    }
    
  }

  EdgeRecord::EdgeRecord(EdgeIter& _edge) : edge(_edge) {

    // TODO: (meshEdit)
    // Compute the combined quadric from the edge endpoints.
    // -> Build the 3x3 linear system whose solution minimizes the quadric error
    //    associated with these two endpoints.
    // -> Use this system to solve for the optimal position, and store it in
    //    EdgeRecord::optimalPoint.
    // -> Also store the cost associated with collapsing this edg in
    //    EdgeRecord::Cost.
    
    Matrix4x4 K = _edge->halfedge()->vertex()->quadric + _edge->halfedge()->twin()->vertex()->quadric;
        	 
    Matrix3x3 A;
        	
    A(0,0) = K(0,0); A(0,1) = K(0,1); A(0,2) = K(0,2);
	A(1,0) = K(1,0); A(1,1) = K(1,1); A(1,2) = K(1,2); 
	A(2,0) = K(2,0); A(2,1) = K(2,1); A(2,2) = K(2,2);  
	
	Vector3D b(-K(0,3), -K(1,3), -K(2,3));
	
	optimalPoint = A.inv() * b;
	
// 	cout<<optimalPoint<<endl;
	
	score = dot(Vector4D(optimalPoint.x, optimalPoint.y, optimalPoint.z, 1.0) , K * Vector4D(optimalPoint.x, optimalPoint.y, optimalPoint.z, 1.0));
	
// 	cout<<score<<endl;
	edge = _edge;

  }

  void MeshResampler::upsample(HalfedgeMesh& mesh)
    // This routine should increase the number of triangles in the mesh using Loop subdivision.
  {

    // TODO: (meshEdit)
    // Compute new positions for all the vertices in the input mesh, using
    // the Loop subdivision rule, and store them in Vertex::newPosition.
    // -> At this point, we also want to mark each vertex as being a vertex of the
    //    original mesh.
    // -> Next, compute the updated vertex positions associated with edges, and
    //    store it in Edge::newPosition.
    // -> Next, we're going to split every edge in the mesh, in any order.  For
    //    future reference, we're also going to store some information about which
    //    subdivided edges come from splitting an edge in the original mesh, and
    //    which edges are new, by setting the flat Edge::isNew. Note that in this
    //    loop, we only want to iterate over edges of the original mesh.
    //    Otherwise, we'll end up splitting edges that we just split (and the
    //    loop will never end!)
    // -> Now flip any new edge that connects an old and new vertex.
    // -> Finally, copy the new vertex positions into final Vertex::position.

    // Each vertex and edge of the original surface can be associated with a vertex in the new (subdivided) surface.
    // Therefore, our strategy for computing the subdivided vertex locations is to *first* compute the new positions
    // using the connectity of the original (coarse) mesh; navigating this mesh will be much easier than navigating
    // the new subdivided (fine) mesh, which has more elements to traverse.  We will then assign vertex positions in
    // the new mesh based on the values we computed for the original mesh.

    // Compute updated positions for all the vertices in the original mesh, using the Loop subdivision rule.

    // Next, compute the updated vertex positions associated with edges.

    // Next, we're going to split every edge in the mesh, in any order.  For future
    // reference, we're also going to store some information about which subdivided
    // edges come from splitting an edge in the original mesh, and which edges are new.
    // In this loop, we only want to iterate over edges of the original mesh---otherwise,
    // we'll end up splitting edges that we just split (and the loop will never end!)

    // Finally, flip any new edge that connects an old and new vertex.

    // Copy the updated vertex positions to the subdivided mesh.
    
//     double u = 0.0f;
    
    for( VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++ )
    {
       v->isNew = false;
       
       double u = (v->degree() == 3) ? 3.f/16.f: 3.f/(8.f * (double)v->degree() );
       v->newPosition = (1.0 - (double)v->degree() * u) * v->position + (u * v->neighborhoodCentroid() * (double)v->degree() );
    }
    
    for( EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++ )
    {
    	e->newPosition = 3.f / 8.f * ( e->halfedge()->vertex()->position + e->halfedge()->twin()->vertex()->position ) + 1.f / 8.f * ( e->halfedge()->next()->twin()->vertex()->position + e->halfedge()->twin()->next()->twin()->vertex()->position ); 
    	e->isNew = false;
    }
    
    
    std::vector<EdgeIter> orig;
    EdgeIter e = mesh.edgesBegin(); 
	while (e != mesh.edgesEnd()) 
	{
		// get the next edge NOW!
		EdgeIter nextEdge = e;
		nextEdge++;

		// now, even if splitting the edge deletes it...
		if (e->isNew == false && e->isBoundary() == false )
		{
			VertexIter v = mesh.splitEdge(e);
			v->position = e->newPosition;
			v->newPosition = e->newPosition;
			v->isNew = true;
			
			v->halfedge()->edge()->isNew = true;
			v->halfedge()->twin()->next()->twin()->next()->edge()->isNew = true;
			
			orig.push_back(v->halfedge()->edge());
			orig.push_back(v->halfedge()->twin()->next()->twin()->next()->edge());
			
			v->halfedge()->twin()->next()->edge()->isNew = true;
			v->halfedge()->next()->next()->edge()->isNew = true;
			
			
		}
		// ...we still have a valid reference to the next edge.
		e = nextEdge;
	}
	
	for( size_t i = 0; i < orig.size(); i++ )
	 {
	 	orig[i]->isNew = false;
	 }
	
	for( EdgeIter e1 = mesh.edgesBegin(); e1 != mesh.edgesEnd(); e1++ )
	{
		if(e1->isNew == true)
		{
			if( (e1->halfedge()->vertex()->isNew == false && e1->halfedge()->twin()->vertex()->isNew == true) || (e1->halfedge()->vertex()->isNew == true && e1->halfedge()->twin()->vertex()->isNew == false))
			{
					mesh.flipEdge(e1);
			}
		}
		
	}

	for( VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++ )
	{
 		v->position = v->newPosition;	
	}
    

  }

  void MeshResampler::downsample(HalfedgeMesh& mesh)
  {

    // TODO: (meshEdit)
    // Compute initial quadrics for each face by simply writing the plane equation
    // for the face in homogeneous coordinates. These quadrics should be stored
    // in Face::quadric
    // -> Compute an initial quadric for each vertex as the sum of the quadrics
    //    associated with the incident faces, storing it in Vertex::quadric
    // -> Build a priority queue of edges according to their quadric error cost,
    //    i.e., by building an EdgeRecord for each edge and sticking it in the
    //    queue.
    // -> Until we reach the target edge budget, collapse the best edge. Remember
    //    to remove from the queue any edge that touches the collapsing edge
    //    BEFORE it gets collapsed, and add back into the queue any edge touching
    //    the collapsed vertex AFTER it's been collapsed. Also remember to assign
    //    a quadric to the collapsed vertex, and to pop the collapsed edge off the
    //    top of the queue.
    
    long int current_num_polygons = 0;
        
    for ( FaceIter f = mesh.facesBegin(); f != mesh.facesEnd(); f++ )
    {
    	f->quadric.zero(0.0);
    	double d = dot(f->normal(), f->halfedge()->vertex()->position);
    	Vector4D v(f->normal().x, f->normal().y, f->normal().z, d);
    	f->quadric = outer(v,v);
    	current_num_polygons++;
//     	cout<<f->quadric<<endl;
    }

    for( VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++ )
    {
    	HalfedgeIter h = v->halfedge();
    	v->quadric.zero(0.0);
    	
    	do{
    		v->quadric = v->quadric + h->face()->quadric;
    		h = h->twin()->next();
    	
    	}while(h != v->halfedge());
    	
//     	cout<<v->quadric<<endl;
    }
        
    MutablePriorityQueue<EdgeRecord> queue;
    
    for( EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++ )
    {
    	EdgeRecord eRecord(e);
    	e->record = eRecord;
    	queue.insert( eRecord );
    }
    
    long int target_num_polygons = current_num_polygons / 4;


	cout<<"current_edges: "<<current_num_polygons<<" "<<"target_edges: "<<target_num_polygons<<endl;

    while(target_num_polygons > 0)
    {
    	printf("Get the best edge\n");
    	//get the best edge record
    	EdgeRecord bestEdge = queue.top();
    	
    	printf("Remove that edge\n");
    	//remove this edge from the queue
    	queue.pop();

		printf("Compute the new quadric\n");
		//compute the new quadric
		Matrix4x4 newK = bestEdge.edge->halfedge()->vertex()->quadric + bestEdge.edge->halfedge()->twin()->vertex()->quadric;
		
		printf("Remove the new edges\n");
		//remove all the edges at the endpoints from the queue
		HalfedgeIter h1 = bestEdge.edge->halfedge()->next();
		HalfedgeIter h2 = bestEdge.edge->halfedge()->twin()->next();
		
		do
		{
			queue.remove(h1->edge()->record);
			h1 = h1->twin()->next();
		}while(h1 != bestEdge.edge->halfedge()->twin());
	
		do
		{
			queue.remove(h2->edge()->record);
			h2 = h2->twin()->next();
		}while(h2 != bestEdge.edge->halfedge());
		
		printf("collapse the new vertex\n");
		//collapse the edge
		VertexIter collapsedVertex = mesh.collapseEdge(bestEdge.edge);
		
		printf("Assign the new quadric to the new vertex\n");
		//assign the new quadric to the collapsed vertex
		collapsedVertex->quadric = newK;
// 		collapsedVertex->position = bestEdge.optimalPoint;
		
		printf("Add all the edges to the queue that touch the new vertex\n");
		//reassign all the edges touching the new vertex and store them again in the queue
		HalfedgeIter h = collapsedVertex->halfedge();
		
// 		printf("Getting the halfedge of the collapsed vertex\n");
		do{
		
		EdgeRecord eRecord(h->edge());
    	h->edge()->record = eRecord;
    	queue.insert( eRecord );
    	h = h->twin()->next();
			
		}while(h != collapsedVertex->halfedge());
		
    	target_num_polygons--;
//     	cout<<target_num_polygons<<endl;
    }
    

  }

  void MeshResampler::resample(HalfedgeMesh& mesh) {

    // TODO: (meshEdit)
    // Compute the mean edge length.
    // Repeat the four main steps for 5 or 6 iterations
    // -> Split edges much longer than the target length (being careful about
    //    how the loop is written!)
    // -> Collapse edges much shorter than the target length.  Here we need to
    //    be EXTRA careful about advancing the loop, because many edges may have
    //    been destroyed by a collapse (which ones?)
    // -> Now flip each edge if it improves vertex degree
    // -> Finally, apply some tangential smoothing to the vertex positions
    
    double meanLength = 0.0;
    int numEdges = 0;
    for ( EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++ )
    {
    	meanLength += e->length();
    	e->isNew = true;
    	numEdges++;
    }
    
    meanLength /= numEdges;
    
    cout<<meanLength<<endl;
    for ( EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++ )
    {  	
    	if(e->length() > 4.0*meanLength/3.0 && e->isNew)
    	{
    		mesh.splitEdge(e);	
    	}
    }
    
    printf("hi\n");
   //  for ( EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++ )
//     {
//     	
//     	if(e->length() < 4.0*meanLength/5.0)
//     	{
//     		mesh.collapseEdge(e);	
//     	}
//     }

	for ( EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++ )
	{
		int a1 = 0, a2 = 0, b1 = 0, b2 = 0;
		
		HalfedgeIter h = e->halfedge()->vertex()->halfedge();
		
		//at edge end points
		do
		{
			a1++;
			h = h->twin()->next();
		}while(h != e->halfedge()->vertex()->halfedge());
		
		h = e->halfedge()->twin()->vertex()->halfedge();
		
		do
		{
			a2++;
			h = h->twin()->next();
		}while(h != e->halfedge()->twin()->vertex()->halfedge());
		
		h = e->halfedge()->next()->twin()->vertex()->halfedge();
		
		//across the edge
		do
		{
			b1++;
			h = h->twin()->next();
		}while(h != e->halfedge()->next()->twin()->vertex()->halfedge());
		
		h = e->halfedge()->twin()->next()->twin()->vertex()->halfedge();
		
		do
		{
			b2++;
			h = h->twin()->next();
		}while(h != e->halfedge()->twin()->next()->twin()->vertex()->halfedge());
		
		int initValence = abs(a1 - 6) + abs(a2 - 6) + abs(b1 - 6) + abs(b2 - 6);
		
		int currValence =  abs(a1 - 1 - 6) + abs(a2 - 1 - 6) + abs(b1 + 1 - 6) + abs(b2 + 1 - 6);
		
		
		if(currValence < initValence)
		{
// 			cout<<"initValence: "<<initValence<<" "<<"currValence: "<<currValence<<endl;
			mesh.flipEdge(e);
		}
		
	}

	for( VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++ )
	{
		v->meanPosition = v->neighborhoodCentroid();
	}
	
	for( VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++ )
	{
		v->position = v->position + (1.0/5.0) * (v->meanPosition - v->position);
	}
    
    

  }

} // namespace CMU462

