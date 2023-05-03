function varargout = trirefine(TR,varargin)

	%% TRIREFINE Uniform triangulation refinement
	%    rTR = trirefine(TR) returns a triangulation that is a uniform
	%    refinement of the input 2-D triangulation TR. Each element of TR is
	%    subdivided into 4 triangles by connecting the midpoints of the edges,
	%    i.e.,
	%
	%             o                              o                  
	%            / \                            / \                     
	%           /   \                          /   \            
	%          /     \         ===>           o-----o
	%         /       \                      / \   / \
	%        /         \                    /   \ /   \
	%       o-----------o                  o-----o-----o
	%     Original triangle          Triangle after refinement
	%
	%    rTRk = trirefine(...,'NumberOfRefinements',k) same as above but 
	%    successively refines the input triangulation TR k times.
	%
	%    [rTR1, rTR2, ..., rTRk] = trirefine(...,'NumberOfRefinements',k) same 
	%    as above but returns the sequence of k refined triangulations. Note: 
	%    The output can also be saved to a cell array using the notation 
	%    [rTR{1:k}] = trirefine(...,'NumberOfRefinements',k).
	%
	%    See also triangulation, triplot, trimesh, trisurf 

	%% Validate input

	if ~isa(TR,'triangulation') || size(TR.ConnectivityList,2)~=3
		error('First input must be a 2-D triangulation object')
	else
		vk = @(x)validateattributes(x,{'numeric'},{'scalar','integer','>=',1});
		ip = inputParser;
		ip.addRequired('TR');
		ip.addParameter('NumberOfRefinements',1,vk);
		ip.parse(TR,varargin{:}); ip.Results;
		k = ip.Results.NumberOfRefinements;
	end

	%% Uniformly refine mesh(es)

	% Initialize the output
	varargout = cell(nargout,1);
	% Loop over the sequence of refinements
	for r = 1:k
		% Retrieve the number of elements and vertices of the triangulation
		nElems = size(TR.ConnectivityList,1);
		nVerts = size(TR.Points,1);
		% Store the (directed "half") edges of the triangulation
		halfEdges{1} = [TR.ConnectivityList(:,1),TR.ConnectivityList(:,2)];
		halfEdges{2} = [TR.ConnectivityList(:,2),TR.ConnectivityList(:,3)];
		halfEdges{3} = [TR.ConnectivityList(:,3),TR.ConnectivityList(:,1)]; 
		% Store the edges of the triangulation
		Edges = TR.edges;
		% Construct the element-to-edge table
		Elem2Edge = zeros(nElems,3);    
		for i = 1:3
			[~,iEdges,iElems]   = intersect([Edges(:,1),Edges(:,2)],...
				halfEdges{i},'rows');
			Elem2Edge(iElems,i) = iEdges;
			[~,iEdges,iElems]   = intersect([Edges(:,2),Edges(:,1)],...
				halfEdges{i},'rows');
			Elem2Edge(iElems,i) = iEdges;
		end
		% Compute and store coordinate of new vertices
		newPoints = (TR.Points(Edges(:,1),:) + TR.Points(Edges(:,2),:))/2;
		% Store new connectivity list for refined elements
		NewConnectivityList1 = [ TR.ConnectivityList(:,1), ...
			nVerts+Elem2Edge(:,1), nVerts+Elem2Edge(:,3) ];
		NewConnectivityList2 = [ nVerts+Elem2Edge(:,1), ...
			TR.ConnectivityList(:,2), nVerts+Elem2Edge(:,2) ];
		NewConnectivityList3 = [ nVerts+Elem2Edge(:,3), ...
			nVerts+Elem2Edge(:,2), TR.ConnectivityList(:,3) ];
		NewConnectivityList4 = nVerts+Elem2Edge(:,1:3);
		% Contruct the refined triangulation
		TR = triangulation([NewConnectivityList1; NewConnectivityList2;...
		  NewConnectivityList3; NewConnectivityList4 ],[TR.Points; newPoints]);
		% Save the output
		varargout{r} = TR;
	end
	if nargout == 1
		varargout{1} = TR;
	end

end