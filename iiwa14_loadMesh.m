function iiwa14_loadMesh()
    %exampleHelperLoadSawyerMesh load sawyer robot mesh for visualization in Simscape Multibody
    %
    %   Set robot's mesh filepaths on the Simscape Multibody visualization blocks.

    %   Copyright 2017-2019 The MathWorks, Inc.
    
    % Define path to Sawyer meshes
    lbrMeshPath = fullfile('./','iiwa_description','meshes','iiwa14','visual');
    
    lbrL0Handle = 'iiwa_paraschos/KUKA LBR iiwa R820 14kg (iiwa14.urdf)/iiwa_link_0/Visual';
    set_param(lbrL0Handle, 'ExtGeomFileName', fullfile(lbrMeshPath, 'link_0.stl')) ;
    
    lbrL1Handle = 'iiwa_paraschos/KUKA LBR iiwa R820 14kg (iiwa14.urdf)/iiwa_link_1/Visual';
    set_param(lbrL1Handle, 'ExtGeomFileName', fullfile(lbrMeshPath, 'link_1.stl')) ;
    
    lbrL2Handle = 'iiwa_paraschos/KUKA LBR iiwa R820 14kg (iiwa14.urdf)/iiwa_link_2/Visual';
    set_param(lbrL2Handle, 'ExtGeomFileName', fullfile(lbrMeshPath, 'link_2.stl')) ;
    
    lbrL3Handle = 'iiwa_paraschos/KUKA LBR iiwa R820 14kg (iiwa14.urdf)/iiwa_link_3/Visual';
    set_param(lbrL3Handle, 'ExtGeomFileName', fullfile(lbrMeshPath, 'link_3.stl')) ;
    
    lbrL4Handle = 'iiwa_paraschos/KUKA LBR iiwa R820 14kg (iiwa14.urdf)/iiwa_link_4/Visual';
    set_param(lbrL4Handle, 'ExtGeomFileName', fullfile(lbrMeshPath, 'link_4.stl')) ;
    
    lbrL5Handle = 'iiwa_paraschos/KUKA LBR iiwa R820 14kg (iiwa14.urdf)/iiwa_link_5/Visual';
    set_param(lbrL5Handle, 'ExtGeomFileName', fullfile(lbrMeshPath, 'link_5.stl')) ;
    
    lbrL6Handle = 'iiwa_paraschos/KUKA LBR iiwa R820 14kg (iiwa14.urdf)/iiwa_link_6/Visual';
    set_param(lbrL6Handle, 'ExtGeomFileName', fullfile(lbrMeshPath, 'link_6.stl')) ;
    
    lbrL7Handle = 'iiwa_paraschos/KUKA LBR iiwa R820 14kg (iiwa14.urdf)/iiwa_link_7/Visual';
    set_param(lbrL7Handle, 'ExtGeomFileName', fullfile(lbrMeshPath, 'link_7.stl')) ;
end