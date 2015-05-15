classdef testReactionTangents < matlab.unittest.TestCase
    % Tests tangent matrices for reaction-diffusion-reaction problems 

    methods (Test)
        function testElementReactionMasterP1(testCase)
            % single element of order 1
            mesh = mkmesh_uniform([0, 1], 1, 1);
            master = mkmaster(mesh);

            jacobians = master.dphi'*mesh.dgnodes;

            % should be the same as the master
            testCase.verifyEqual(jacobians, ones(size(jacobians)), 'Abstol', 1e-10);

            b = 1;
            A = element_reaction_tangent(master, jacobians, b);
            testCase.verifyEqual(A, (b/6)*[2 1; 1 2], 'Abstol', 1e-10);

            b = 24.2;
            A = element_reaction_tangent(master, jacobians, b);
            testCase.verifyEqual(A, (b/6)*[2 1; 1 2], 'Abstol', 1e-10);
        end

        function testElementReactionP1(testCase)
            % single element of order 1
            mesh = mkmesh_uniform([0, 0.5], 1, 1);
            master = mkmaster(mesh);

            jacobians = master.dphi'*mesh.dgnodes;

            % should be the half the size of master
            testCase.verifyEqual(jacobians, 0.5*ones(size(jacobians)), 'Abstol', 1e-10);

            % reaction tangents should also halve
            b = 1;
            A = element_reaction_tangent(master, jacobians, b);
            testCase.verifyEqual(A, 0.5*(b/6)*[2 1; 1 2], 'Abstol', 1e-10);

            b = 192.6;
            A = element_reaction_tangent(master, jacobians, b);
            testCase.verifyEqual(A, 0.5*(b/6)*[2 1; 1 2], 'Abstol', 1e-10);
        end
    end
end
