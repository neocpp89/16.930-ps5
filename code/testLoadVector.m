classdef testLoadVector < matlab.unittest.TestCase
    % Tests tangent matrices for convection-diffusion-reaction problems 

    methods (Test)
        function testLoadVectorMasterP1(testCase)
            % single element of order 1
            mesh = mkmesh_uniform([0, 1], 1, 1);
            master = mkmaster(mesh);

            jacobians = master.dphi'*mesh.dgnodes;

            % should be the same as the master
            testCase.verifyEqual(jacobians, ones(size(jacobians)), 'Abstol', 1e-10);

            fn = @(x) 1.2*x+20;

            fg = fn(master.gp);
            f = element_load_vector(master, jacobians, fg);

            % from exact integration
            testCase.verifyEqual(f, [10.2; 10.4], 'Abstol', 1e-10);
        end

        function testLoadVectorP1(testCase)
            % single element of order 1
            mesh = mkmesh_uniform([0, 0.5], 1, 1);
            master = mkmaster(mesh);

            jacobians = master.dphi'*mesh.dgnodes;

            % should be the half the size of master
            testCase.verifyEqual(jacobians, 0.5*ones(size(jacobians)), 'Abstol', 1e-10);

            fn = @(x) 1.2*x+20;

            gp = master.phi'*mesh.dgnodes;

            fg = fn(gp);
            f = element_load_vector(master, jacobians, fg);

            % from exact integration
            testCase.verifyEqual(f, [5.05; 5.1], 'Abstol', 1e-10);
        end
    end
end
