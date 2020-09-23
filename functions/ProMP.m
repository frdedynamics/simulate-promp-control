classdef ProMP %< handle
    %Class definition for a ProMP Datatype
    %   > in progress <
    % help and reference in "Representing Structured Data with Classes"
    % matlab documentation
    % web(fullfile(docroot, 'matlab/matlab_oop/example-representing-structured-data.html'))
    
    % Indexing: promp.property(NrOfDemo)
    
    properties
        Data MPData %data needs to be a MPData object
        Phase
        Label = ['MyProMP @ ' datestr(datetime('now'))]
        BasisFcns
        BfBlockDiag_full
        BfBlockDiag_tvec
        BasisFcns_dot
        Bf_dotBlockDiag_full
        Bf_dotBlockDiag_tvec
        BasisFcns_dotdot
        Bf_dotdotBlockDiag_full
        Bf_dotdotBlockDiag_tvec   
        Bf_and_Bf_dotBlockDiag_tvec
        Bf_dot_and_Bf_dotdotBlockDiag_tvec
        NrBasisFcns = 10 % default 10 %set the SetAccess as immutable such that it can only be set in the constructor
        NrDOF
        NrDemos
        NrTimesteps
        Trainer ProMPTrainer
        Weights MPWeights
        Mean
        Variance
        Variance_MLE
        CovY
        %TODO: make properties private as necessary
    end
        
    methods
        %% ProMP Constructor
        function promp = ProMP(NrBasisFcns, NrDOF)
            %CONSTRUCTOR of the class to easily generate instances of it
            if nargin == 2
                promp.NrBasisFcns = NrBasisFcns;
                promp.NrDOF = NrDOF;
            else
                error('Pleas specify Number of Basis Functions and Number of DOF.')
            end
        end
            
        %% User Methods
        function promp = trainProMP(promp)
            [trainer, trainedProMP] = train(promp.Trainer, promp);            
            promp = trainedProMP;
            promp.Trainer = trainer;
        end
        
        function promp = addData(promp, newData, trainNow)
            %ADD DATA: add new demonstration to the ProMP
            %TODO: Align demos in time by Dynamic time Warping if neccessary

            %% make sure new data is a MPData object:
            newData = MPData(newData.Time, newData.q); 
            
            %% Resample all demos if desired
            resample = false; %true            
            if resample
                resampled_NrTimesteps = 250;
                resampled_duration = 1; %seconds
                q_adjustedLength = interp1(linspace(0,1,length(newData.q)), newData.q, linspace(0,1,resampled_NrTimesteps)');
                Time_adjustedLength = linspace(0, resampled_duration, resampled_NrTimesteps)';
                newData = MPData(Time_adjustedLength, q_adjustedLength);
            end
            
            %% Handling the first demo:
            if numel(promp.Data) == 0 
                promp.Data(1) = newData; %add first dataset
                promp.NrDemos = size(promp.Data,2);
                promp = adjustOutputPhase(promp, 0.01, 1); % generate Phase signal according to default values dt=0.01 seconds and duration 1 seconds

                promp.CovY = eye(promp.NrDOF); % define dimensions of CovY matrix and initialize with identity matrix
                promp.Variance = eye(promp.NrDOF * promp.NrBasisFcns);
                promp.Mean = zeros(promp.NrBasisFcns, promp.NrDOF);
            else
            %% Handling all subsequent demos:       
                promp.Data = [promp.Data newData];  %add new data set    

                % Update NrOfDemos property
                promp.NrDemos = size(promp.Data,2);
            end
            
            %% Handle training method input
            if nargin < 3 
                trainNow = 'trainBat'; %set default value
            end
            
            if trainNow == 'trainNow'
                promp = trainProMP(promp);
            elseif trainNow == 'trainBat'
                %Batch training has to be called manually after adding all
                %the data
            end
            
        end
        
        function promp = adjustOutputPhase(promp, dt, duration)
            % Generate phase signal with desired output dt and duration

            NrTimesteps = (1/dt) * duration;
            promp.Phase = linspace(0, duration, NrTimesteps);
            
            
            %Update Nr of Timesteps
            promp.NrTimesteps = NrTimesteps;
            
            %Update the BasisFunctions to match the desired output phase
            promp = generateOutputBasisFcns(promp, promp.Phase);
        end
                   
        function qdot_prior = generateQdotPrior(promp)
            % Idea: abs(diff(bf) * promp.Weights(1).Weights) % where bf is a
            % basisfunction matrix with as many timesteps as there are
            % weights and as many basis functions as the promp has
            lambda = 10e-12; %ridge regression regularization

            for i=1:promp.NrDOF
                firstDemoWeights_acc(:,i) = (promp.BasisFcns'*promp.BasisFcns + lambda*eye(size(promp.BasisFcns'*promp.BasisFcns,1)))\promp.BasisFcns'*[0;0; diff(promp.Data(end).q(:,i),2)];
            end  
                        
            for i=1:promp.NrDOF
                firstDemoWeights(:,i) = (promp.BasisFcns'*promp.BasisFcns + lambda*eye(size(promp.BasisFcns'*promp.BasisFcns,1)))\promp.BasisFcns'*promp.Data(end).q(:,i);
            end            
            
            N = promp.NrBasisFcns; %   N = number of basis functions
            T = promp.NrBasisFcns; %+ 1; 
            h = 0.5*(T/(N-1))^2; % h = 0.5*(c_{i+1} - c{i})^2, also found in Rückert Probabilistic Machine Learning
            
            weight_size_Bf = gaussianBasisFcns(N, T, h); % scale the basis function matrix down from NrOfTimesteps length to NrofWeights length
            weight_size_Bf_dot = gradient(weight_size_Bf', 1/T)';%[zeros(1,N); diff(weight_size_Bf)];
            weight_size_Bf_dot_dot = gradient(weight_size_Bf_dot', 1/T)'; %[zeros(1,N); diff(weight_size_Bf_dot)];
            %qdot_weight_size = abs(weight_size_Bf * firstDemoWeights_acc);
            qdot_weight_size = abs(weight_size_Bf_dot * firstDemoWeights);
            %qdot_weight_size = abs(weight_size_Bf_dot_dot * firstDemoWeights); %*0.005
            qdot_prior = blockDiagEye(promp.NrDOF, promp.NrBasisFcns) .* (qdot_weight_size(:)*qdot_weight_size(:)'); 
            %qdot_prior = eye(promp.NrDOF * promp.NrBasisFcns) .* qdot_weight_size(:); 
            %qdot_prior = qdot_prior ./ max(max(qdot_prior)); % normalizing
        end
        
        function promp = conditionToWP(promp, WP, CovWP, phase_t)
            %CONDITION TO WAYPOINT: Condition the trajectory distribution
            %to a certain waypoint and a corresponding uncertainty in the
            %waypoint. According to Paraschos18
            % Create Block diagonal BasisFunction matrix for all DoFs at
            % time step t.
            t = round(phase_t * (length(promp.Phase)-1) +1); %convert from phase to time step
            NrOfBlocks = promp.NrDOF;
            BfBlockDiag_t = kron(eye(NrOfBlocks),promp.BasisFcns(t,:))'; %the transpose of the kronecker product is strange, dimensions should be correct without it according to paraschos18 but they aren't
            L = promp.Variance * BfBlockDiag_t / (CovWP + BfBlockDiag_t' * promp.Variance * BfBlockDiag_t);
            %%%%%%L = promp.Variance * promp.BasisFcns(t,:) / (CovWP + promp.BasisFcns(t,:)' * promp.Variance * promp.BasisFcns(t,:));
            mean_cond = promp.Mean(:) + L*(WP - BfBlockDiag_t' * promp.Mean(:));
            cov_cond = promp.Variance - L * BfBlockDiag_t' * promp.Variance;
            promp.Mean = reshape(mean_cond, promp.NrBasisFcns, promp.NrDOF);
            promp.Variance = cov_cond;
        end
        
        function plotEncodedData(promp, StartIdx, EndIdx, MultipleFigures)
            %PLOT ENCODED DATA: plot the encoded trajectories to check the
            % quality of the regression. Requires the Index of the first
            % and the last Demo to plot. MultipleFigures Parameter is
            % binary and switches between plotting in separate figures(=1) or a
            % single one (=0).
            if MultipleFigures == 1
                % Plot Data of each Demonstration in seperate figure
                for i=StartIdx:EndIdx
                    figure;
                    plot(promp.Phase, promp.BasisFcns * promp.Weights(i).Weights);
                end
                legend('Joint 1', 'Joint 2', 'Joint 3', 'Joint 4', 'Joint 5', 'Joint 6', 'Joint 7'); title('Joint Positions');                
            elseif MultipleFigures == 0
                % Plot Data of all demos in single figure. Good to see
                % variance
                figure; hold on
                for i=StartIdx:EndIdx                    
                    plot(promp.Phase, promp.BasisFcns * promp.Weights(i).Weights);
                end
                % Plot the mean trajectory:
                %hold on;
                mean_trajectory = promp.BasisFcns * promp.Mean;
                plot(promp.Phase, mean_trajectory,'LineWidth',2); 
                legend('Joint 1', 'Joint 2', 'Joint 3', 'Joint 4', 'Joint 5', 'Joint 6', 'Joint 7'); title('Joint Positions');
            else
                disp('Error: Specify if you want multiple (1) or a single (0) figure');
            end
        end
        
        function plotMeanVarTraj(promp, subplotDOFs)
            %PLOT THE ENCODED MEAN AND STANDARD DEVIATION TRAJECTORIES
            
            % Handle input arguments
            if nargin < 2
                subplotDOFs = 'SinglePlot';
            end
            
            % Generate Mean trajectories
            % States/Positions
            mean_trajectory = promp.BasisFcns * promp.Mean;
                                   
            % margingal distribution of p(y) in lin gauss system y=Ax+b.
            % See Murphy Equation 4.126            
            % Generate Variance trajectories
            cov_trajectory = zeros(promp.NrDOF, promp.NrDOF, promp.NrTimesteps);
            for t=1:promp.NrTimesteps
                cov_trajectory(:,:,t) = promp.CovY + promp.BfBlockDiag_tvec(:,:,t)*promp.Variance*promp.BfBlockDiag_tvec(:,:,t)'; % promp.CovY + 
            end
            
            % Plot
            MeanVarFig = figure; %to use this figure again: figure(MeanVarFig)            
            cmap = colormap(lines);
            if strcmpi(subplotDOFs, 'subplotDOFs')  
                set(MeanVarFig,'units','normalized','outerposition',[0 0 0.28 1]);
                for k = 1:promp.NrDOF
                    figure(MeanVarFig)
                    %q_var(:,k) = cov_trajectory(k,k,:); %for saving the data
                    subplot(promp.NrDOF, 1, k)
                    s=shadedErrorBar(promp.Phase, mean_trajectory(:,k), 2*sqrt(cov_trajectory(k,k,:)));
                    s.patch.FaceColor = cmap(k,:);
                    s.mainLine.Color = cmap(k,:);
                    s.mainLine.LineWidth = 1;
                    s.edge(1).LineStyle = 'none';
                    s.edge(2).LineStyle = 'none';
                    legend(['DOF ' num2str(k)]);
                end 
                sgtitle('Mean-Variance Trajectories');
                
            elseif strcmpi(subplotDOFs, 'SinglePlot')
                figure(MeanVarFig); hold on;
                for k = 1:promp.NrDOF                    
                    %q_var(:,k) = cov_trajectory(k,k,:); %for saving the data                    
                    s=shadedErrorBar(promp.Phase, mean_trajectory(:,k), 2*sqrt(cov_trajectory(k,k,:)));
                    s.patch.FaceColor = cmap(k,:);
                    s.mainLine.Color = cmap(k,:);
                    s.mainLine.LineWidth = 1;
                    s.edge(1).LineStyle = 'none';
                    s.edge(2).LineStyle = 'none';
                    legendEntries{k} = ['DOF ' num2str(k)];
                end
                legend(legendEntries); title('Mean-Variance Trajectories');
            end
        end
        
        function plotMeanVarTrajDot(promp, subplotDOFs)
            %PLOT THE ENCODED MEAN AND STANDARD DEVIATION TRAJECTORIES
            
            % Handle input arguments
            if nargin < 2
                subplotDOFs = 'SinglePlot';
            end
            
            % Generate Mean trajectories            
            % Time derivatives of states / Velocities
            mean_trajectory_dot = promp.BasisFcns_dot * promp.Mean;
                     
            % margingal distribution of p(y) in lin gauss system y=Ax+b.
            % See Murphy Equation 4.126            
            % Generate Variance trajectories
            for t=1:promp.NrTimesteps
                cov_trajectory_dot(:,:,t) = promp.Bf_dotBlockDiag_tvec(:,:,t)*promp.Variance*promp.Bf_dotBlockDiag_tvec(:,:,t)';
            end
            
            % Plot
            MeanVarFig = figure; %to use this figure again: figure(MeanVarFig)            
            cmap = colormap(lines);
            if strcmpi(subplotDOFs, 'subplotDOFs')  
                set(MeanVarFig,'units','normalized','outerposition',[0 0 0.28 1]);
                for k = 1:promp.NrDOF
                    figure(MeanVarFig)
                    %q_var(:,k) = cov_trajectory_dot(k,k,:); %for saving the data
                    subplot(promp.NrDOF, 1, k)
                    s=shadedErrorBar(promp.Phase, mean_trajectory_dot(:,k), 2*sqrt(cov_trajectory_dot(k,k,:)));
                    s.patch.FaceColor = cmap(k,:);
                    s.mainLine.Color = cmap(k,:);
                    s.mainLine.LineWidth = 1;
                    s.edge(1).LineStyle = 'none';
                    s.edge(2).LineStyle = 'none';
                    legend(['DOF ' num2str(k)]);
                end 
                sgtitle('Mean-Variance Trajectories 1. Derivative');
                
            elseif strcmpi(subplotDOFs, 'SinglePlot')
                figure(MeanVarFig); hold on;
                for k = 1:promp.NrDOF                    
                    %q_var(:,k) = cov_trajectory_dot(k,k,:); %for saving the data                    
                    s=shadedErrorBar(promp.Phase, mean_trajectory_dot(:,k), 2*sqrt(cov_trajectory_dot(k,k,:)));
                    s.patch.FaceColor = cmap(k,:);
                    s.mainLine.Color = cmap(k,:);
                    s.mainLine.LineWidth = 1;
                    s.edge(1).LineStyle = 'none';
                    s.edge(2).LineStyle = 'none';
                    legendEntries{k} = ['DOF ' num2str(k)];
                end
                legend(legendEntries); title('Mean-Variance Trajectories 1. Derivative');
            end
        end        

        function plotMeanVarTrajDotDot(promp, subplotDOFs)
            %PLOT THE ENCODED MEAN AND STANDARD DEVIATION TRAJECTORIES
            
            % Handle input arguments
            if nargin < 2
                subplotDOFs = 'SinglePlot';
            end
            
            % Generate Mean trajectories            
            % Time derivatives of states / Velocities
            mean_trajectory_dotdot = promp.BasisFcns_dotdot * promp.Mean;
                     
            % margingal distribution of p(y) in lin gauss system y=Ax+b.
            % See Murphy Equation 4.126            
            % Generate Variance trajectories
            for t=1:promp.NrTimesteps
                cov_trajectory_dotdot(:,:,t) = promp.Bf_dotdotBlockDiag_tvec(:,:,t)*promp.Variance*promp.Bf_dotdotBlockDiag_tvec(:,:,t)';
            end
            
            % Plot
            MeanVarFig = figure; %to use this figure again: figure(MeanVarFig)            
            cmap = colormap(lines);
            if strcmpi(subplotDOFs, 'subplotDOFs')  
                set(MeanVarFig,'units','normalized','outerposition',[0 0 0.28 1]);
                for k = 1:promp.NrDOF
                    figure(MeanVarFig)
                    %q_var(:,k) = cov_trajectory_dot(k,k,:); %for saving the data
                    subplot(promp.NrDOF, 1, k)
                    s=shadedErrorBar(promp.Phase, mean_trajectory_dotdot(:,k), 2*sqrt(cov_trajectory_dotdot(k,k,:)));
                    s.patch.FaceColor = cmap(k,:);
                    s.mainLine.Color = cmap(k,:);
                    s.mainLine.LineWidth = 1;
                    s.edge(1).LineStyle = 'none';
                    s.edge(2).LineStyle = 'none';
                    legend(['DOF ' num2str(k)]);
                end 
                sgtitle('Mean-Variance Trajectories 2. Derivative');
                
            elseif strcmpi(subplotDOFs, 'SinglePlot')
                figure(MeanVarFig); hold on;
                for k = 1:promp.NrDOF                    
                    %q_var(:,k) = cov_trajectory_dot(k,k,:); %for saving the data                    
                    s=shadedErrorBar(promp.Phase, mean_trajectory_dotdot(:,k), 2*sqrt(cov_trajectory_dotdot(k,k,:)));
                    s.patch.FaceColor = cmap(k,:);
                    s.mainLine.Color = cmap(k,:);
                    s.mainLine.LineWidth = 1;
                    s.edge(1).LineStyle = 'none';
                    s.edge(2).LineStyle = 'none';
                    legendEntries{k} = ['DOF ' num2str(k)];
                end
                legend(legendEntries); title('Mean-Variance Trajectories 2. Derivative');
            end
        end        

        function plot3MeanVarTraj(promp)
            mean_trajectory = promp.BasisFcns * promp.Mean;
            for t=1:promp.NrTimesteps
                cov_trajectory(:,:,t) = promp.CovY + promp.BfBlockDiag_tvec(:,:,t)*promp.Variance*promp.BfBlockDiag_tvec(:,:,t)';
            end
            
            for i = 1:3
                q_var(:,i) = cov_trajectory(i,i,:);
            end
            figure;
            grid on; 
            hold on;
            view(3);
            %axis([0 1 0 1 0 1]); 
            axis equal
            plot3(mean_trajectory(:,1),mean_trajectory(:,2),mean_trajectory(:,3), 'LineWidth', 1.5);
            tubeplot_r3(mean_trajectory(:,1),mean_trajectory(:,2),mean_trajectory(:,3),2*sqrt(q_var),2*sqrt(mean(q_var,2)),30);
           
            % Fancy lighting
            camlight
            lighting gouraud
            %material dull
            shading interp
            
            xlabel('X'),ylabel('Y'),zlabel('Z');
            view(-136.289,36.7318);
            %keep in mind that the renderer settings change how the
            %horizontal lines look. painters looks nicer but is slower
            %TODO : add orientation
        end
        
        function plot3TaskSpaceMeanVarTraj(promp)
            robot = loadrobot('kukaIiwa14', 'DataFormat', 'row'); % 'kukaIiwa14' % 'universalUR5' % 'frankaEmikaPanda'
            
            % Joint-Space mean trajectory
            mean_trajectory = promp.BasisFcns * promp.Mean;            
            
            homog_transform = zeros(4,4,promp.NrTimesteps);
            mean_trajectory_translation_ts = zeros(promp.NrTimesteps, 3);
            mean_trajectory_quat_ts = zeros(promp.NrTimesteps, 4);
            
            for t = 1:promp.NrTimesteps
                % Forward kinematics: get homogeneous transformations from base to flange for all joint configurations (joint-space mean trajectory)
                homog_transform(:,:,t) = getTransform(robot, mean_trajectory(t,:), 'iiwa_link_ee');
                
                % extract translations out of the homog. transformation matrices
                mean_trajectory_translation_ts(t,:) = homog_transform(1:3,4,t);
                
                % convert rotations in homog. trans. matrices to quaternions (for plotTransforms() ) 
                mean_trajectory_quat_ts(t,:) = tform2quat(homog_transform(:,:,t));
            end
           
            % Linearize forward kinematics by first order tayler expansion --> jacobian
            for t=1:promp.NrTimesteps
                jacobian(:,:,t) = geometricJacobian(robot, mean_trajectory(t,:), 'iiwa_link_ee');
                cov_trajectory(:,:,t) =  promp.BfBlockDiag_tvec(:,:,t)*promp.Variance*promp.BfBlockDiag_tvec(:,:,t)';%+ promp.CovY ;
                cov_trajectory_ts(:,:,t) = jacobian(:,:,t)*cov_trajectory(:,:,t)*jacobian(:,:,t)';
            end
            
            % Get the Variance for translation in X,Y,Z from the task space
            % covariance matrices
            for i = 4:6
                q_var(:,i-3) = cov_trajectory_ts(i,i,:);
            end
            
            % Plot
            set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
            view(3);
            %axis([ -0.4 0.4 -0.7 -0.3 0.2 0.8])
            
            %1) plot coordinate frames to visualize orientations
            NrFramesToPlot = promp.NrTimesteps / ((8e-10)*promp.NrTimesteps^3 + 18);
            frame_idx = round(linspace(1, promp.NrTimesteps, NrFramesToPlot));
            plotTransforms(mean_trajectory_translation_ts(frame_idx,:), mean_trajectory_quat_ts(frame_idx,:), 'FrameSize',0.1);
            hold on
            grid on
            
            %2) plot mean in 3D
            plot3(mean_trajectory_translation_ts(:,1),mean_trajectory_translation_ts(:,2),mean_trajectory_translation_ts(:,3), 'k', 'LineWidth', 1.5);
            
            %3) plot Variance as tube
            NrRowsToPlot = promp.NrTimesteps / ((1e-10)*promp.NrTimesteps^3 + 1);
            row_idx = round(linspace(1,promp.NrTimesteps,NrRowsToPlot));
            tubeplot_r3(mean_trajectory_translation_ts(row_idx,1),mean_trajectory_translation_ts(row_idx,2),mean_trajectory_translation_ts(row_idx,3),2*sqrt(q_var(row_idx,:)),2*sqrt(mean(q_var(row_idx,:),2)),30, [1 1 1]);            
            axis image
            hold off
            
            % Fancy lighting
            camlight
            lighting gouraud
            %material dull
            shading interp
            
            %orbit camera once:
%             for i=1:36*2
%                 camorbit(5,0)                
%                 drawnow
%                 pause(0.01)
%             end
            
            xlabel('X'),ylabel('Y'),zlabel('Z');
            view(-136.289,36.7318);
        end
        
        function plotMeanTrajRobot(promp)
            %% Corke Toolbox Approach
%                 DH PARAMETERS FOR THE ROBOT iiwa 14 R 820
%                 alpha=[pi/2 -pi/2 -pi/2 pi/2 pi/2 -pi/2 0]';
%                 d=[0.36 0.0 0.42 0.0 0.4 0.0 0.126]';            
%                 a=[0 0 0 0 0 0 0]';
%                 theta=[pi 0 0 0 0 0 0]';            
%                 dh = [theta d a alpha];
%                 %% Generate SerialLink object from DH params (Corke Toolbox)
%                 rob = SerialLink(dh);
%                 view(-145.5,30.0)
%                 rob.plot([0 0 0 0 0 0 0])
%                 %rob.teach
%                 %% Compute Mean trajectory
%                 mean_trajectory = promp.BasisFcns * promp.Mean;
%                 delay = 1/promp.NrTimesteps;
%                 for t=1:25:promp.NrTimesteps
%                     rob.plot(mean_trajectory(t,:))
%                     %pause(delay)
%                 end

            %% Matlab Robotics Systems toolbox
            figure;
            plot3TaskSpaceMeanVarTraj(promp)
            hold on
            robot = importrobot('iiwa14.urdf');
            robot.DataFormat = 'row';
            show(robot, 'PreservePlot',false,'Frames', 'off');
            axis([-1 1 -1 1 0 1.5]);
            view_angles = [61.26,24.98]; % [45.18,23.95]
            view(view_angles);            
            camlight %new camlight that at position of the current view
          
            mean_trajectory = promp.BasisFcns * promp.Mean;
            %pause(0.1)
            increment = round((2e-10)*promp.NrTimesteps^3 + 10);
            for t=1:increment:promp.NrTimesteps                
                show(robot, mean_trajectory(t,:),'PreservePlot',false,'Frames', 'off');
                lighting gouraud
                drawnow
                %axis image
                axis([-1 1 -1 1 0 1.5]);
                view(view_angles)
                %pause(0.0005)
                t
                
                if t > promp.NrTimesteps - increment
                    show(robot, mean_trajectory(end,:),'PreservePlot',false,'Frames', 'off');
                    lighting gouraud
                    t = promp.NrTimesteps
                    %axis image
                    axis([-1 1 -1 1 0 1.5]);
                    view(view_angles)
                    break
                end
                
            end
        end
        
        function [q_mean, q_var, q_covar] = computeMeanVarTraj(promp)
            q_mean = promp.BasisFcns * promp.Mean;
            q_mean_dot = promp.BasisFcns_dot * promp.Mean;
            q_mean_dotdot = promp.BasisFcns_dotdot * promp.Mean;
                                        
            % Generate Variance trajectories
            q_covar = zeros(promp.NrDOF, promp.NrDOF, promp.NrTimesteps);
            q_var = zeros(promp.NrTimesteps, promp.NrDOF);
            for t=1:promp.NrTimesteps
                q_covar(:,:,t) = promp.CovY + promp.BfBlockDiag_tvec(:,:,t)*promp.Variance*promp.BfBlockDiag_tvec(:,:,t)'; % promp.CovY +
                q_var(t,:) = diag(q_covar(:,:,t));
            end          
            
        end
          
        function [k, K, epsilon_u] = computeControlGains(promp)
            A = blockDiagMat(promp.NrDOF, [0 1; 0 0]);
            B = blockDiagMat(promp.NrDOF, [0; 1]);
            c = zeros(promp.NrDOF*2, 1);
                       
            Sigma_w = promp.Variance;
            mu_w = promp.Mean(:);            
            
            for t = 1:promp.NrTimesteps-1
                % time step dependent stuff:
                
                % BasisFcns
                Psi_t = promp.Bf_and_Bf_dotBlockDiag_tvec(:,:,t);
                Psi_tplusdt = promp.Bf_and_Bf_dotBlockDiag_tvec(:,:,t+1);
                
                % variance for timestep t
                Sigma_t = Psi_t * Sigma_w * Psi_t';
                Sigma_tplusdt = Psi_tplusdt * Sigma_w * Psi_tplusdt';
                
                %% Control Noise
                C_t = Psi_t * Sigma_w * Psi_tplusdt';
                Sigma_s(:,:,t) = (Sigma_tplusdt - (C_t' * (Sigma_t \ C_t)))./ (promp.Phase(t+1)-promp.Phase(t));
                Sigma_s(:,:,t) = (Sigma_s(:,:,t) + Sigma_s(:,:,t)') * 0.5 +  eye(size(Sigma_s(:,:,t)))*1e-12;
                Sigma_u(:,:,t) = B \ (Sigma_s(:,:,t)) / B'; 
                Sigma_u(:,:,t) = (Sigma_u(:,:,t) + Sigma_u(:,:,t)') * 0.5;
                epsilon_u(:,t) = mvnrnd(zeros(promp.NrDOF,1), Sigma_u(:,:,t)./(promp.Phase(t+1)-promp.Phase(t)));
            end
            
            for t = 1:promp.NrTimesteps-1
                % time step dependent stuff:
                
                % BasisFcns
                Psi_dot_t = promp.Bf_dot_and_Bf_dotdotBlockDiag_tvec(:,:,t);
                Psi_t = promp.Bf_and_Bf_dotBlockDiag_tvec(:,:,t);
                
                % variance for timestep t
                Sigma_t = Psi_t * Sigma_w * Psi_t';
                
                %% Feedback gain matrix
                K(:,:,t) = B \ (Psi_dot_t * Sigma_w * Psi_t' - A * Sigma_t - 0.5 * Sigma_s(:,:,t)) / Sigma_t;
                
                %     disp('K: '); disp(K);
                
                %% Feedforward gains
                k(:,t) = B \ (Psi_dot_t * mu_w - (A + B*K(:,:,t)) * Psi_t * mu_w - c);
                %     disp('k: '); disp(k);
            end
        end
        
        function sampledDemos = sampleDemosFromProMP(promp, NrSamples)
            % Samples a Demonstration from the Trajectory Distribution
            % encoded a ProMP
            if nargin < 2
                NrSamples = 1;
            end
            
            for n = 1:NrSamples
                % Sample a weight vector first
                sampledWeights = mvnrnd(promp.Mean(:), promp.Variance)';

                % Sample states from the distribution parametrized by the
                % sampled weight vector and CovY
                sampledTraj = zeros(promp.NrTimesteps, promp.NrDOF);
                for t=1:promp.NrTimesteps
                    sampledTraj(t,:) = mvnrnd(promp.BfBlockDiag_tvec(:,:,t)*sampledWeights, promp.CovY);
                end    
                
                sampledDemo = MPData(promp.Phase, sampledTraj);
                sampledDemos(n) = sampledDemo;
            
                % Plot
                % figure;plot(sampledDemo.q);
            end
        end
                
    end
    
    %% Private Methods
    methods ( Access = private ) 
        
        function promp = generateOutputBasisFcns(promp, output_time_vector)
            % GENERATE OUTPUT BASIS FUNCTIONS: generates the basis function matrices used for generating mean-variance trajectories,
            % plotting, sampling and computing controller gains.       
                        
            BfStruct = generateBasisFcns(promp.NrBasisFcns, promp.NrDOF, output_time_vector);
            
            promp.BasisFcns = BfStruct.BasisFcns;
            promp.BfBlockDiag_full = BfStruct.BfBlockDiag_full;
            promp.BfBlockDiag_tvec = BfStruct.BfBlockDiag_tvec;
            BasisFcns_unnormalized = BfStruct.BasisFcns_unnormalized;
            z_dot = BfStruct.z_dot;
            
            %% Derivatives of BasisFcns
            promp.BasisFcns_dot = gradient(promp.BasisFcns', promp.Phase, ones(length(promp.NrBasisFcns),1))'; %* z_dot; %*z_dot not needed here ?! since time not phase dependent?!
            promp.BasisFcns_dotdot = gradient(promp.BasisFcns_dot', promp.Phase, ones(length(promp.NrBasisFcns),1))';
            
            promp.Bf_dotBlockDiag_full = blockDiagMat(promp.NrDOF, promp.BasisFcns_dot);
            promp.Bf_dotdotBlockDiag_full = blockDiagMat(promp.NrDOF, promp.BasisFcns_dotdot);
            
            promp.Bf_dotBlockDiag_tvec = zeros(promp.NrDOF, promp.NrDOF*promp.NrBasisFcns, length(output_time_vector));
            promp.Bf_dotdotBlockDiag_tvec = zeros(promp.NrDOF, promp.NrDOF*promp.NrBasisFcns, length(output_time_vector));
            for t = 1:length(output_time_vector)
                promp.Bf_dotBlockDiag_tvec(:,:,t) = blockDiagMat(promp.NrDOF, promp.BasisFcns_dot(t,:));
                promp.Bf_dotdotBlockDiag_tvec(:,:,t) = blockDiagMat(promp.NrDOF, promp.BasisFcns_dotdot(t,:));                
                
                % Basis function matrix containing Bf and Bf_dot. This for
                % is needed to compute the control gains of the stochstic
                % feedback controller in paraschos
                promp.Bf_and_Bf_dotBlockDiag_tvec(:,:,t) = blockDiagMat(promp.NrDOF, [promp.BasisFcns(t,:); promp.BasisFcns_dot(t,:)]);
                promp.Bf_dot_and_Bf_dotdotBlockDiag_tvec(:,:,t) = blockDiagMat(promp.NrDOF, [promp.BasisFcns_dot(t,:); promp.BasisFcns_dotdot(t,:)]);
            end
            
        end
        
    end   
end




