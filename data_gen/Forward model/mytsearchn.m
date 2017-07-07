function [ind, int_func] = mytsearchn(mesh,coord)

% Determines in which mesh element a given point falls. Similar to matlab tsearchn function 

if mesh.dimension == 2
	[N,junk] = size(coord);
	ind = NaN(N,1); int_func = NaN(N,3);
	for i = 1:N
		dist = (mesh.nodes(:,1:2)-repmat(coord(i,1:2),size(mesh.nodes(:,1:2),1),1)).^2;
		dist = sqrt(dist(:,1) + dist(:,2));
		[snodes,temp_ind] = sort(dist);
		clear dist
		j=1; true = 0;
		while true == 0 && j <= 10
			[r,c] = find(mesh.elements == temp_ind(j));
			foo = mesh.elements(r,:);
			[n,m] = size(foo);
			k = 1; true = 0;
			while true == 0 & k <= n
				P = mesh.nodes(foo(k,1),1:2); Q = mesh.nodes(foo(k,2),1:2); R = mesh.nodes(foo(k,3),1:2);
				true=in_tri(coord(i,1:2),P,Q,R);
				if true == 1
					ind(i,1) = r(k);
					A = [P(1) Q(1) R(1); P(2) Q(2) R(2); 1 1 1];
					b = [coord(i,1); coord(i,2); 1];
					int_func(i,:) = (A\b)';
				elseif true == 0
					k = k+1;
				end
			end
			j=j+1;
		end
	end

elseif mesh.dimension == 3

	[N,junk] = size(coord);
	ind = NaN(N,1); int_func = NaN(N,4);
	for i = 1:N
			dist = (mesh.nodes(:,1:3)-repmat(coord(i,1:3),size(mesh.nodes(:,1:3),1),1)).^2;
			dist = sqrt(dist(:,1) + dist(:,2) + dist(:,3));
		[snodes,temp_ind] = sort(dist);
		clear dist
		j=1; true = 0;
		while true == 0 && j <= 10
			[r,c] = find(mesh.elements == temp_ind(j));
			foo = mesh.elements(r,:);
			[n,m] = size(foo);
			k = 1; true = 0;
			while true == 0 && k <= n
				P = mesh.nodes(foo(k,1),1:3); 
				Q = mesh.nodes(foo(k,2),1:3); 
				R = mesh.nodes(foo(k,3),1:3);
				S = mesh.nodes(foo(k,4),1:3);
				true=in_tetra(coord(i,1:3),P,Q,R,S);
				if true == 1
					ind(i,1) = r(k);
					A = [P(1) Q(1) R(1) S(1); P(2) Q(2) R(2) S(2); P(3) Q(3) R(3) S(3); 1 1 1 1];
					b = [coord(i,1); coord(i,2); coord(i,3); 1];
					int_func(i,:) = (A\b)';
				elseif true == 0
					k = k+1;
				end
			end
			j=j+1;
		end
	end

end


%%%%%%%%%%%%%% sub function space %%%%%%%%%%%%%%%%%%
function True=in_tri(P,P1,P2,P3)

%in_tri is used to check if a point P is inside the triangle P1P2P3 or not. 


Area_P1P2P3 = 1/2. *abs(det([P1(1) P1(2) 1;P2(1) P2(2) 1;P3(1) P3(2) 1]));

if abs(1/2. *((abs(det([P(1) P(2) 1;P1(1) P1(2) 1;P2(1) P2(2) 1])))...
        + (abs(det([P(1) P(2) 1;P2(1) P2(2) 1;P3(1) P3(2) 1])))...
        + (abs(det([P(1) P(2) 1;P3(1) P3(2) 1;P1(1) P1(2) 1])))...
        )-Area_P1P2P3)/Area_P1P2P3 < 10^-6
    True = 1;
else
    True = 0;
end

function True=in_tetra(P,P1,P2,P3,P4)

%in_tetra is used to check if a point P is inside the tetrahedron P1P2P3P4 or not. 

V0 = 1/6.*abs(det([P1(1) P1(2) P1(3) 1;P2(1) P2(2) P2(3) 1;P3(1) P3(2) P3(3) 1;P4(1) P4(2) P4(3) 1]));

if abs(1/6.*((abs(det([P(1) P(2) P(3) 1;P2(1) P2(2) P2(3) 1;P3(1) P3(2) P3(3) 1;P4(1) P4(2) P4(3) 1])))...
        + (abs(det([P1(1) P1(2) P1(3) 1;P(1) P(2) P(3) 1;P3(1) P3(2) P3(3) 1;P4(1) P4(2) P4(3) 1])))...
        + (abs(det([P1(1) P1(2) P1(3) 1;P2(1) P2(2) P2(3) 1;P(1) P(2) P(3) 1;P4(1) P4(2) P4(3) 1])))...
        + (abs(det([P1(1) P1(2) P1(3) 1;P2(1) P2(2) P2(3) 1;P3(1) P3(2) P3(3) 1;P(1) P(2) P(3) 1])))...
        ) - V0)/V0 < 10^-6
    True = 1;
else
    True = 0;
end

