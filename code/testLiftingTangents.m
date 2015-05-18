classdef testLiftingTangents < matlab.unittest.TestCase
    % Tests tangent matrices for convection-diffusion-reaction problems 

    methods (Test)
        function testFaceLiftingMasterP1(testCase)
            % single element of order 1
            mesh = mkmesh_uniform([0, 1], 1, 1);
            master = mkmaster(mesh);

            jacobians = master.dphi'*mesh.dgnodes;

            % should be the same as the master
            testCase.verifyEqual(jacobians, ones(size(jacobians)), 'Abstol', 1e-10);

            nu = 1;
            A = face_lifting_tangent(master, master, jacobians, jacobians, nu);
            testCase.verifyEqual(A, 4*nu*[1, -1; -1 1], 'Abstol', 1e-10);

            nu = 0.3;
            A = face_lifting_tangent(master, master, jacobians, jacobians, nu);
            testCase.verifyEqual(A, 4*nu*[1, -1; -1 1], 'Abstol', 1e-10);
        end

        function testFaceLiftingP1(testCase)
            % single element of order 1
            mesh = mkmesh_uniform([0, 0.5], 1, 1);
            master = mkmaster(mesh);

            jacobians = master.dphi'*mesh.dgnodes;

            % should be the half the size of master
            testCase.verifyEqual(jacobians, 0.5*ones(size(jacobians)), 'Abstol', 1e-10);

            % tangent components double 
            nu = 1;
            A = face_lifting_tangent(master, master, jacobians, jacobians, nu);
            testCase.verifyEqual(A, 8*nu*[1, -1; -1 1], 'Abstol', 1e-10);

            nu = 0.3;
            A = face_lifting_tangent(master, master, jacobians, jacobians, nu);
            testCase.verifyEqual(A, 8*nu*[1, -1; -1 1], 'Abstol', 1e-10);
        end
    end
end
