classdef testConvectionTangents < matlab.unittest.TestCase
    % Tests tangent matrices for convection-diffusion-reaction problems 

    methods (Test)
        function testElementConvectionMasterP1(testCase)
            % single element of order 1
            mesh = mkmesh_uniform([0, 1], 1, 1);
            master = mkmaster(mesh);

            jacobians = master.dphi'*mesh.dgnodes;

            % should be the same as the master
            testCase.verifyEqual(jacobians, ones(size(jacobians)), 'Abstol', 1e-10);

            c = 1;
            A = element_convection_tangent(master, jacobians, c);
            testCase.verifyEqual(A, 0.5*c*[1 1; -1 -1], 'Abstol', 1e-10);

            c = 0.5;
            A = element_convection_tangent(master, jacobians, c);
            testCase.verifyEqual(A, 0.5*c*[1 1; -1 -1], 'Abstol', 1e-10);

            c = 15.2134;
            A = element_convection_tangent(master, jacobians, c);
            testCase.verifyEqual(A, 0.5*c*[1 1; -1 -1], 'Abstol', 1e-10);
        end

        function testElementConvectionP1(testCase)
            % single element of order 1
            mesh = mkmesh_uniform([0, 0.5], 1, 1);
            master = mkmaster(mesh);

            jacobians = master.dphi'*mesh.dgnodes;

            % should be the half the size of master
            testCase.verifyEqual(jacobians, 0.5*ones(size(jacobians)), 'Abstol', 1e-10);

            % but convection tangents do not change
            c = 1;
            A = element_convection_tangent(master, jacobians, c);
            testCase.verifyEqual(A, 0.5*c*[1 1; -1 -1], 'Abstol', 1e-10);

            c = 0.5;
            A = element_convection_tangent(master, jacobians, c);
            testCase.verifyEqual(A, 0.5*c*[1 1; -1 -1], 'Abstol', 1e-10);

            c = 15.2134;
            A = element_convection_tangent(master, jacobians, c);
            testCase.verifyEqual(A, 0.5*c*[1 1; -1 -1], 'Abstol', 1e-10);
        end
    end
end
