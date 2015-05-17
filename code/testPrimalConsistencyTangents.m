classdef testPrimalConsistencyTangents < matlab.unittest.TestCase
    % Tests tangent matrices for convection-diffusion-reaction problems 

    methods (Test)
        function testFacePrimalConsistencyMasterP1(testCase)
            % single element of order 1
            mesh = mkmesh_uniform([0, 1], 1, 1);
            master = mkmaster(mesh);

            jacobians = master.dphi'*mesh.dgnodes;

            % should be the same as the master
            testCase.verifyEqual(jacobians, ones(size(jacobians)), 'Abstol', 1e-10);

            nu = 1;
            [AL, AR] = face_primal_consistency_tangent(master, master, jacobians, jacobians, nu);
            testCase.verifyEqual(AL, 0.5*nu*[1, -1; -1 1], 'Abstol', 1e-10);
            testCase.verifyEqual(AR, 0.5*nu*[1, -1; -1 1], 'Abstol', 1e-10);

            nu = 0.3;
            [AL, AR] = face_primal_consistency_tangent(master, master, jacobians, jacobians, nu);
            testCase.verifyEqual(AL, 0.5*nu*[1, -1; -1 1], 'Abstol', 1e-10);
            testCase.verifyEqual(AR, 0.5*nu*[1, -1; -1 1], 'Abstol', 1e-10);
        end

        function testFacePrimalConsistencyP1(testCase)
            % single element of order 1
            mesh = mkmesh_uniform([0, 0.5], 1, 1);
            master = mkmaster(mesh);

            jacobians = master.dphi'*mesh.dgnodes;

            % should be the half the size of master
            testCase.verifyEqual(jacobians, 0.5*ones(size(jacobians)), 'Abstol', 1e-10);

            % tangent components double 
            nu = 1;
            [AL, AR] = face_primal_consistency_tangent(master, master, jacobians, jacobians, nu);
            testCase.verifyEqual(AL, nu*[1, -1; -1 1], 'Abstol', 1e-10);
            testCase.verifyEqual(AR, nu*[1, -1; -1 1], 'Abstol', 1e-10);

            nu = 0.3;
            [AL, AR] = face_primal_consistency_tangent(master, master, jacobians, jacobians, nu);
            testCase.verifyEqual(AL, nu*[1, -1; -1 1], 'Abstol', 1e-10);
            testCase.verifyEqual(AR, nu*[1, -1; -1 1], 'Abstol', 1e-10);
        end
    end
end
