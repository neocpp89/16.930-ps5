classdef testDiffusionTangents < matlab.unittest.TestCase
    % Tests tangent matrices for convection-diffusion-reaction problems 

    methods (Test)
        function testElementDiffusionMasterP1(testCase)
            % single element of order 1
            mesh = mkmesh_uniform([0, 1], 1, 1);
            master = mkmaster(mesh);

            jacobians = master.dphi'*mesh.dgnodes;

            % should be the same as the master
            testCase.verifyEqual(jacobians, ones(size(jacobians)), 'Abstol', 1e-10);

            nu = 1;
            A = element_diffusion_tangent(master, jacobians, nu);
            testCase.verifyEqual(A, nu*[1 -1; -1 1], 'Abstol', 1e-10);

            nu = 125.1;
            A = element_diffusion_tangent(master, jacobians, nu);
            testCase.verifyEqual(A, nu*[1 -1; -1 1], 'Abstol', 1e-10);
        end

        function testElementDiffusionP1(testCase)
            % single element of order 1
            mesh = mkmesh_uniform([0, 0.5], 1, 1);
            master = mkmaster(mesh);

            jacobians = master.dphi'*mesh.dgnodes;

            % should be the half the size of master
            testCase.verifyEqual(jacobians, 0.5*ones(size(jacobians)), 'Abstol', 1e-10);

            % tangents components should double
            nu = 1;
            A = element_diffusion_tangent(master, jacobians, nu);
            testCase.verifyEqual(A, 2*nu*[1 -1; -1 1], 'Abstol', 1e-10);

            nu = 104;
            A = element_diffusion_tangent(master, jacobians, nu);
            testCase.verifyEqual(A, 2*nu*[1 -1; -1 1], 'Abstol', 1e-10);
        end
    end
end
