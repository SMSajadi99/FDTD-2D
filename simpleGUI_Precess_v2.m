function simpleGUI_Precess_v2
    % Create a figure window for the GUI
    fig = uifigure('Position', [500, 300, 400, 200], 'Name', 'Welcome PHOTORA-v0.2');

    % Add a welcome label to the GUI
    uilabel(fig, ...
        'Position', [100, 120, 200, 50], ...
        'Text', 'Welcome to the Program!', ...
        'FontSize', 16, ...
        'HorizontalAlignment', 'center');

    % Add a button to start the program
    uibutton(fig, ...
        'Position', [150, 70, 100, 30], ...
        'Text', 'Start', ...
        'ButtonPushedFcn', @(btn, event) showWaveguideSelection(fig));

    % Add your name at the bottom of the GUI
    uilabel(fig, ...
        'Position', [100, 20, 200, 20], ...
        'Text', 'Created by Seyed Mohammad Sajadi', ...
        'FontSize', 10, ...
        'HorizontalAlignment', 'center');
end

function showWaveguideSelection(prevFig)
    % Close the previous GUI
    close(prevFig);
    
    % Create a new figure for waveguide selection
    fig = uifigure('Name', 'Waveguide Selection', 'Position', [250 100 900 450]);

    % Create a button group without a title and border
    bg = uibuttongroup(fig, 'Position', [50 350 800 50], 'Title', '', 'BorderType', 'none');

    % Radio buttons inside the group with larger font size
    rb1 = uiradiobutton(bg, 'Text', 'RIB Waveguide', 'Position', [60 10 150 25], 'FontSize', 14);
    rb2 = uiradiobutton(bg, 'Text', 'Planer Waveguide', 'Position', [320 10 160 25], 'FontSize', 14);
    rb3 = uiradiobutton(bg, 'Text', 'Square Waveguide', 'Position', [590 10 180 25], 'FontSize', 14);

    % Axes for displaying images (shifted left)
    ax1 = uiaxes(fig, 'Position', [50 120 250 200]);
    ax2 = uiaxes(fig, 'Position', [325 120 250 200]);
    ax3 = uiaxes(fig, 'Position', [600 120 250 200]);

    % Confirm button
    btn = uibutton(fig, 'Text', 'Confirm', 'FontSize', 14, ...
                   'Position', [400 50 100 50], ...
                   'ButtonPushedFcn', @(btn,event) confirmSelection(bg));

    % Add separator lines
    line1 = uilabel(fig, 'Position', [305 100 2 300], 'BackgroundColor', 'black');
    line2 = uilabel(fig, 'Position', [580 100 2 300], 'BackgroundColor', 'black');

    % Load images with transparency effect
    loadTransparentImage(ax1, './images/rib.png');
    loadTransparentImage(ax2, './images/planer.png');
    loadTransparentImage(ax3, './images/square.png');
end

function loadTransparentImage(ax, imgPath)
    if exist(imgPath, 'file')
        img = imread(imgPath); % Read image

        % Convert to grayscale to detect white background
        grayImg = rgb2gray(img);
        mask = grayImg < 240; % Threshold to detect non-white regions
        
        % Convert mask to alpha data (Transparency)
        alpha = double(mask);

        % Display image with transparency
        imshow(img, 'Parent', ax);
        set(ax.Children, 'AlphaData', alpha);
    else
        uialert(ax.Parent, ['Image file not found: ', imgPath], 'Error');
    end
end

% Define the function to handle waveguide selection confirmation
function confirmSelection(bg)
    % Get the selected waveguide type
    selectedButton = bg.SelectedObject;
    if isempty(selectedButton)
        uialert(bg.Parent, 'Please select a waveguide type.', 'Selection Required');
    else
        waveguideType = selectedButton.Text;
        disp(['You selected: ', waveguideType]);
        if strcmp(waveguideType, 'RIB Waveguide')
            showRIBInputs();
        elseif strcmp(waveguideType, 'Planer Waveguide')
            showPlanerInputs();
        else
            showRectangularInputs();
        end
    end
end
  
% Define the function to show inputs for RIB Wave-guide
function showRIBInputs()
    % Create a new figure for RIB inputs
    fig = uifigure('Position', [80, 50, 1200, 550], 'Name', 'RIB Wave-guide Inputs');

    % Add an image (enlarged and positioned)
    axImage = uiaxes(fig, 'Position', [450, 140, 750, 400]); % Larger position
    imshow('./images/rib_prosess.png', 'Parent', axImage, 'Border', 'tight'); % Remove border
    axImage.XColor = 'none';
    axImage.YColor = 'none';

    % Labels and input fields for refractive indices
    uilabel(fig, 'Position', [20, 500, 200, 20], 'Text', 'Refractive Index of Air:');
    airIndex = uieditfield(fig, 'numeric', 'Position', [280, 500, 100, 20], 'Value', 1);

    uilabel(fig, 'Position', [20, 470, 200, 20], 'Text', 'Refractive Index of Waveguide:');
    waveguideIndex = uieditfield(fig, 'numeric', 'Position', [280, 470, 100, 20], 'Value', 3.3);

    uilabel(fig, 'Position', [20, 440, 200, 20], 'Text', 'Refractive Index of Substrate:');
    substrateIndex = uieditfield(fig, 'numeric', 'Position', [280, 440, 100, 20], 'Value', 1.5);

    % Labels and input fields for overall dimensions
    uilabel(fig, 'Position', [20, 400, 250, 20], 'Text', 'Overall Dimensions (X Start, X End):');
    overallXStart = uieditfield(fig, 'numeric', 'Position', [280, 400, 50, 20], 'Value', 0);
    overallXEnd = uieditfield(fig, 'numeric', 'Position', [340, 400, 50, 20], 'Value', 10);

    uilabel(fig, 'Position', [20, 370, 250, 20], 'Text', 'Overall Dimensions (Y Start, Y End):');
    overallYStart = uieditfield(fig, 'numeric', 'Position', [280, 370, 50, 20], 'Value', 0);
    overallYEnd = uieditfield(fig, 'numeric', 'Position', [340, 370, 50, 20], 'Value', 10);

    % Labels and input fields for substrate dimensions
    uilabel(fig, 'Position', [20, 330, 250, 20], 'Text', 'Substrate Dimensions (X Start, X End):');
    substrateXStart = uieditfield(fig, 'numeric', 'Position', [280, 330, 50, 20], 'Value', 0);
    substrateXEnd = uieditfield(fig, 'numeric', 'Position', [340, 330, 50, 20], 'Value', 10);

    uilabel(fig, 'Position', [20, 300, 250, 20], 'Text', 'Substrate Dimensions (Y Start, Y End):');
    substrateYStart = uieditfield(fig, 'numeric', 'Position', [280, 300, 50, 20], 'Value', 0);
    substrateYEnd = uieditfield(fig, 'numeric', 'Position', [340, 300, 50, 20], 'Value', 6);

    % Labels and input fields for waveguide dimensions
    uilabel(fig, 'Position', [20, 260, 250, 20], 'Text', 'Waveguide Dimensions (X Start, X End):');
    waveguideXStart = uieditfield(fig, 'numeric', 'Position', [280, 260, 50, 20], 'Value', 0);
    waveguideXEnd = uieditfield(fig, 'numeric', 'Position', [340, 260, 50, 20], 'Value', 10);

    uilabel(fig, 'Position', [20, 230, 250, 20], 'Text', 'Waveguide Dimensions (Y Start, Y End):');
    waveguideYStart = uieditfield(fig, 'numeric', 'Position', [280, 230, 50, 20], 'Value', 6);
    waveguideYEnd = uieditfield(fig, 'numeric', 'Position', [340, 230, 50, 20], 'Value', 9);

    % Labels and input fields for Waveguide RIB Dimensions (X)
    uilabel(fig, 'Position', [20, 190, 250, 20], 'Text', 'Waveguide RIB Dimensions (X Start, X End):');
    ribXStart = uieditfield(fig, 'numeric', 'Position', [280, 190, 50, 20], 'Value', 3.5);
    ribXEnd = uieditfield(fig, 'numeric', 'Position', [340, 190, 50, 20], 'Value', 6.5);
    
    % Labels and input fields for Waveguide RIB Dimensions (Y)
    uilabel(fig, 'Position', [20, 160, 250, 20], 'Text', 'Waveguide RIB Dimensions (Y Start, Y End):');
    ribYStart = uieditfield(fig, 'numeric', 'Position', [280, 160, 50, 20], 'Value', 9);
    ribYEnd = uieditfield(fig, 'numeric', 'Position', [340, 160, 50, 20], 'Value', 10);

    % Labels and input fields for mesh points
    uilabel(fig, 'Position', [20, 120, 200, 20], 'Text', 'Number of Mesh Points (X, Y):');
    meshPointsX = uieditfield(fig, 'numeric', 'Position', [280, 120, 50, 20], 'Value', 49);
    meshPointsY = uieditfield(fig, 'numeric', 'Position', [340, 120, 50, 20], 'Value', 49);

    % New fields for Mode and Lambda
    uilabel(fig, 'Position', [20, 80, 200, 20], 'Text', 'Mode:');
    modeField = uieditfield(fig, 'numeric', 'Position', [280, 80, 100, 20], 'Value', 1);

    uilabel(fig, 'Position', [20, 50, 200, 20], 'Text', 'Lambda (Wavelength):');
    lambdaField = uieditfield(fig, 'numeric', 'Position', [280, 50, 100, 20], 'Value', 1.55);
    
    % Create a styled label for the final answer box
    uilabel(fig, 'Position', [800, 130, 200, 30], 'Text', 'Confinement Factor', ...
        'FontSize', 14, 'FontWeight', 'bold', 'FontColor', [0 0 0.5]); % Dark blue color
    
    % Create a box for the final answer at the bottom right
    finalAnswerBox = uitextarea(fig, 'Position', [760, 100, 200, 30], 'Editable', 'on', ...
        'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', [0.9 0.9 0.9]); % Light gray background

    % Add a confirm button (Note: pass finalAnswerBox to the callback)
    uibutton(fig, 'Position', [200, 10, 100, 30], 'Text', 'Submit', ...
        'ButtonPushedFcn', @(btn, event) submitRIBInputs(airIndex, waveguideIndex, substrateIndex, ...
        overallXStart, overallXEnd, overallYStart, overallYEnd, ...
        substrateXStart, substrateXEnd, substrateYStart, substrateYEnd, ...
        waveguideXStart, waveguideXEnd, waveguideYStart, waveguideYEnd, ...
        ribXStart, ribXEnd, ribYStart, ribYEnd, ...
        meshPointsX, meshPointsY, modeField, lambdaField, finalAnswerBox));
end


function submitRIBInputs(airIndex, waveguideIndex, substrateIndex, overallXStart, overallXEnd, overallYStart, overallYEnd, ...
    substrateXStart, substrateXEnd, substrateYStart, substrateYEnd, waveguideXStart, waveguideXEnd, waveguideYStart, waveguideYEnd, ...
    ribXStart, ribXEnd, ribYStart, ribYEnd, meshPointsX, meshPointsY, modeField, lambdaField, finalAnswerBox)

    % Gather input values into a structure (if needed)
    inputs = struct( ...
        'AirIndex', airIndex.Value, ...
        'WaveguideIndex', waveguideIndex.Value, ...
        'SubstrateIndex', substrateIndex.Value, ...
        'OverallX', [overallXStart.Value, overallXEnd.Value], ...
        'OverallY', [overallYStart.Value, overallYEnd.Value], ...
        'SubstrateX', [substrateXStart.Value, substrateXEnd.Value], ...
        'SubstrateY', [substrateYStart.Value, substrateYEnd.Value], ...
        'WaveguideX', [waveguideXStart.Value, waveguideXEnd.Value], ...
        'WaveguideY', [waveguideYStart.Value, waveguideYEnd.Value], ...
        'WaveguideRIBX', [ribXStart.Value, ribXEnd.Value], ...
        'WaveguideRIBY', [ribYStart.Value, ribYEnd.Value], ...
        'MeshPoints', [meshPointsX.Value, meshPointsY.Value], ...
        'ModeField', modeField.Value, ...
        'LambdaField', lambdaField.Value ...
    );

    disp('RIB Wave-guide Inputs:');
    disp(inputs);
    
    % ---------------------- Your Computational Code Starts Here ---------------------- %
    clc;
    % (Assuming you do not want to close all figures during computation)
    % close all;
    
    % Define domain sizes based on overall dimensions
    Lx = overallXEnd.Value - overallXStart.Value;
    Ly = overallYEnd.Value - overallYStart.Value;
    
    % Number of mesh points
    Nx = meshPointsX.Value;
    Ny = meshPointsY.Value;
    
    dx = Lx / Nx;
    dy = Ly / Ny;
    
    % Generate grid coordinates
    x = linspace(0, Lx, Nx + 1); % x-coordinates
    y = linspace(0, Ly, Ny + 1); % y-coordinates
    
    % Initialize the refractive index matrix with air (n = 1.0)
    n_matrix = ones(Nx - 1, Ny - 1);
    
    % Define substrate region
    substrate_x_min = substrateXStart.Value;
    substrate_x_max = substrateXEnd.Value;
    substrate_y_min = substrateYStart.Value;
    substrate_y_max = substrateYEnd.Value;
    
    % Set substrate refractive index in its region
    for i = 1:(Nx - 1)
        for j = 1:(Ny - 1)
            if x(i) >= substrate_x_min && x(i) <= substrate_x_max && ...
               y(j) >= substrate_y_min && y(j) <= substrate_y_max
                n_matrix(i, j) = substrateIndex.Value;
            end
        end
    end
    
    % Process for RIB waveguide
    waveguide_type = 'RIB';
    
    if strcmpi(waveguide_type, 'RIB')
        % Base rectangle of the waveguide
        base_x_min = waveguideXStart.Value;
        base_x_max = waveguideXEnd.Value;
        base_y_min = waveguideYStart.Value;
        base_y_max = waveguideYEnd.Value;
    
        for i = 1:(Nx - 1)
            for j = 1:(Ny - 1)
                if x(i) >= base_x_min && x(i) <= base_x_max && ...
                   y(j) >= base_y_min && y(j) <= base_y_max
                    n_matrix(i, j) = waveguideIndex.Value;
                end
            end
        end
    
        % Rib on top of the base
        rib_x_min = ribXStart.Value;
        rib_x_max = ribXEnd.Value;
        rib_y_min = ribYStart.Value;
        rib_y_max = ribYEnd.Value;
    
        for i = 1:(Nx - 1)
            for j = 1:(Ny - 1)
                if x(i) >= rib_x_min && x(i) <= rib_x_max && ...
                   y(j) >= rib_y_min && y(j) <= rib_y_max
                    n_matrix(i, j) = waveguideIndex.Value;
                end
            end
        end
    end
    
    % Visualize the refractive index distribution (optional)
    figure;
    imagesc(x, y, n_matrix');
    colorbar;
    title('n(x, y) Distribution');
    xlabel('x-axis');
    ylabel('y-axis');
    set(gca, 'YDir', 'normal');
    
    % Define wavelength and compute k0
    lambda = lambdaField.Value; % Wavelength
    k0 = 2 * pi / lambda;
    disp(['k0 = ', num2str(k0)]);
    
    % Setup the sparse matrix for eigenvalue problem
    deltaX = dx;
    deltaY = dy;
    N_total = (Nx - 1) * (Ny - 1);
    A = sparse(N_total, N_total);
    
    % Populate matrix A (finite difference approximation)
    for i = 1:(Nx - 1)
        for j = 1:(Ny - 1)
            index = (j - 1) * (Nx - 1) + i;
            n_value = n_matrix(i, j);
            A(index, index) = -2 / (deltaX^2) - 2 / (deltaY^2) + (k0^2) * (n_value^2);
            if j > 1
                A(index, index - (Nx - 1)) = 1 / (deltaX^2);
            end
            if j < Ny - 1
                A(index, index + (Nx - 1)) = 1 / (deltaX^2);
            end
            if i > 1
                A(index, index - 1) = 1 / (deltaY^2);
            end
            if i < Nx - 1
                A(index, index + 1) = 1 / (deltaY^2);
            end
        end
    end
    
    % Compute eigenvalues and eigenvectors
    disp('Computing eigenvalues and eigenvectors...');
    [eigenvectors, eigenvalues_matrix] = eigs(A, N_total);
    eigenvalues = diag(eigenvalues_matrix);
    
    % Use modeField.Value as the index for selecting the eigenvalue/eigenvector
    index = modeField.Value;
    while index < 1 || index > N_total
        disp('Invalid index. Please select a valid index.');
        index = input(['Select an eigenvalue index (1 to ', num2str(N_total), '): ']);
    end
    
    selected_eigenvalue = eigenvalues(index);
    selected_eigenvector = eigenvectors(:, index);
    
    % Plot the selected eigenvector as a mode shape
    figure;
    mode_shape = reshape(selected_eigenvector, [Nx - 1, Ny - 1]);
    imagesc(abs(mode_shape)');
    colorbar;
    title(['Eigenvector (Mode Shape) for Eigenvalue ', num2str(selected_eigenvalue)]);
    xlabel('x-axis');
    ylabel('y-axis');
    set(gca, 'YDir', 'normal');
    
    % -------------------- Confinement Factor Calculation -------------------- %
    % Compute the total energy (sum of |E|^2 over the entire domain)
    total_energy = sum(sum(abs(mode_shape).^2));
    
    % Create a mask for the waveguide region using the provided waveguide dimensions.
    % Note: x and y are defined for grid points; use x(1:end-1) and y(1:end-1) to match mode_shape.
    waveguide_mask = (x(1:end-1) >= waveguideXStart.Value & x(1:end-1) <= waveguideXEnd.Value)' & ...
                     (y(1:end-1) >= waveguideYStart.Value & y(1:end-1) <= waveguideYEnd.Value);
    
    % Compute the energy confined in the waveguide region
    waveguide_energy = sum(sum(abs(mode_shape(waveguide_mask)).^2));
    
    % Calculate the Confinement Factor
    confinement_factor = waveguide_energy / total_energy;
    
    % Display the Confinement Factor in the console
    disp(['Confinement Factor: ', num2str(confinement_factor)]);
    
    % Show the result in the final answer box in the GUI
    finalAnswerBox.Value = num2str(confinement_factor);
end

% Define the function to show inputs for Planer Wave-guide
function showPlanerInputs()
    % Create a new figure for Planer inputs
    fig = uifigure('Position', [20, 100, 1200, 450], 'Name', 'Planer Waveguide Inputs');

    % Add an image (enlarged and positioned)
    axImage = uiaxes(fig, 'Position', [450, 40, 750, 400]); % Larger position
    imshow('./images/planer_prosess.png', 'Parent', axImage, 'Border', 'tight'); % Remove border
    axImage.XColor = 'none';
    axImage.YColor = 'none';


    % Labels and input fields for refractive indices
    uilabel(fig, 'Position', [20, 400, 200, 20], 'Text', 'Refractive Index of Air:');
    airIndex = uieditfield(fig, 'numeric', 'Position', [280, 400, 100, 20], 'Value', 1);

    uilabel(fig, 'Position', [20, 370, 200, 20], 'Text', 'Refractive Index of Waveguide:');
    waveguideIndex = uieditfield(fig, 'numeric', 'Position', [280, 370, 100, 20], 'Value', 3.3);

    uilabel(fig, 'Position', [20, 340, 200, 20], 'Text', 'Refractive Index of Substrate:');
    substrateIndex = uieditfield(fig, 'numeric', 'Position', [280, 340, 100, 20], 'Value', 1.5);

    % Labels and input fields for overall dimensions
    uilabel(fig, 'Position', [20, 300, 200, 20], 'Text', 'Overall Length:');
    overallLength = uieditfield(fig, 'numeric', 'Position', [280, 300, 100, 20], 'Value', 10);

    uilabel(fig, 'Position', [20, 270, 200, 20], 'Text', 'Overall Width:');
    overallWidth = uieditfield(fig, 'numeric', 'Position', [280, 270, 100, 20], 'Value', 10);

    % Labels and input fields for thicknesses
    uilabel(fig, 'Position', [20, 230, 200, 20], 'Text', 'Thickness of Waveguide:');
    waveguideThickness = uieditfield(fig, 'numeric', 'Position', [280, 230, 100, 20], 'Value', 2);

    uilabel(fig, 'Position', [20, 200, 200, 20], 'Text', 'Thickness of Substrate:');
    substrateThickness = uieditfield(fig, 'numeric', 'Position', [280, 200, 100, 20], 'Value', 4);

    % Labels and input fields for mesh points
    uilabel(fig, 'Position', [20, 160, 200, 20], 'Text', 'Number of Mesh Points (X, Y):');
    meshPointsX = uieditfield(fig, 'numeric', 'Position', [280, 160, 50, 20], 'Value', 49);
    meshPointsY = uieditfield(fig, 'numeric', 'Position', [340, 160, 50, 20], 'Value', 49);

    % Labels and input fields for mode and lambda
    uilabel(fig, 'Position', [20, 120, 200, 20], 'Text', 'Mode:');
    modeField = uieditfield(fig, 'numeric', 'Position', [280, 120, 100, 20], 'Value', 1);

    uilabel(fig, 'Position', [20, 90, 200, 20], 'Text', 'Lambda (Wavelength):');
    lambdaField = uieditfield(fig, 'numeric', 'Position', [280, 90, 100, 20], 'Value', 1.55);
    
    % Create a styled label for the final answer box
    uilabel(fig, 'Position', [800, 60, 200, 30], 'Text', 'Confinement Factor', ...
        'FontSize', 14, 'FontWeight', 'bold', 'FontColor', [0 0 0.5]); % Dark blue color
    
    % Create a box for the final answer at the bottom right
    finalAnswerBox = uitextarea(fig, 'Position', [760, 30, 200, 30], 'Editable', 'on', ...
        'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', [0.9 0.9 0.9]); % Light gray background

    % Add a confirm button
    uibutton(fig, 'Position', [200, 30, 100, 30], 'Text', 'Submit', ...
    'ButtonPushedFcn', @(btn, event) submitPlanerInputs(airIndex, waveguideIndex, substrateIndex, ...
    overallLength, overallWidth, waveguideThickness, substrateThickness, ...
    meshPointsX, meshPointsY, modeField, lambdaField, finalAnswerBox));
end

% Function to handle Planer inputs submission
function submitPlanerInputs(airIndex, waveguideIndex, substrateIndex, overallLength, overallWidth, waveguideThickness,...
    substrateThickness, meshPointsX, meshPointsY, modeField, lambdaField, finalAnswerBox)
    % Gather input values
    inputs = struct( ...
        'AirIndex', airIndex.Value, ...
        'WaveguideIndex', waveguideIndex.Value, ...
        'SubstrateIndex', substrateIndex.Value, ...
        'OverallLength', overallLength.Value, ...
        'OverallWidth', overallWidth.Value, ...
        'WaveguideThickness', waveguideThickness.Value, ...
        'SubstrateThickness', substrateThickness.Value, ...
        'MeshPoints', [meshPointsX.Value, meshPointsY.Value], ...
        'ModeField', modeField.Value, ...
        'LambdaField', lambdaField.Value ...
    );

    disp('Planer Wave-guide Inputs:');
    disp(inputs);
    % Proceed with further processing or computation
    
    % --- Your Computational Code Starts Here ---
    clc;
    close all;
    
    % Input values for domain and mesh size
    Lx = overallLength.Value;
    Ly = overallWidth.Value;
    
    % Calculate the number of mesh points
    Nx = meshPointsX.Value;
    Ny = meshPointsY.Value;
    
    dx = Lx/Nx;
    dy = Ly/Ny;
    
    % Generate grid coordinates
    x = linspace(0, Lx, Nx + 1); % x-coordinates
    y = linspace(0, Ly, Ny + 1); % y-coordinates
    
    % Initialize the n(x, y) matrix with air (n = 1.0)
    n_matrix = ones(Nx - 1, (Ny - 1));
    
    % Input substrate region
    substrate_x_min = 0;
    substrate_x_max = overallLength.Value;
    substrate_y_min = 0;
    substrate_y_max = substrateThickness.Value;
    
    % Assign refractive index for the substrate (n = 1.5)
    for i = 1:(Nx - 1)
        for j = 1:(Ny - 1)
            if x(i) >= substrate_x_min && x(i) <= substrate_x_max && ...
               y(j) >= substrate_y_min && y(j) <= substrate_y_max
                n_matrix(i, j) = substrateIndex.Value;
            end
        end
    end
    
    % Input waveguide parameters
    waveguide_type = 'Plannar';
    
    if strcmpi(waveguide_type, 'Plannar')
        % Input for base rectangle of waveguide
        base_x_min = 0;
        base_x_max = overallLength.Value;
        base_y_min = substrateThickness.Value;
        base_y_max = waveguideThickness.Value + substrateThickness.Value;
    
        % Assign refractive index for the base rectangle (n = 3.3)
        for i = 1:(Nx - 1)
            for j = 1:(Ny - 1)
                if x(i) >= base_x_min && x(i) <= base_x_max && ...
                   y(j) >= base_y_min && y(j) <= base_y_max
                    n_matrix(i, j) = waveguideIndex.Value;
                end
            end
        end
    end
    
    % Display the resulting n(x, y) matrix
    % disp('n(x, y) matrix:');
    % disp(n_matrix);
    
    % Visualize the grid
    figure;
    imagesc(x, y, n_matrix');
    colorbar;
    title('n(x, y) Distribution');
    xlabel('x-axis');
    ylabel('y-axis');
    
    % Reverse the y-axis direction
    set(gca, 'YDir', 'normal');

    lambda = lambdaField.Value; % Wavelength
    k0 = 2 * pi / lambda; % k0 calculation
    % disp(['k0 = ', num2str(k0)]);
    
    % Mesh size (from previous step)
    deltaX = dx;
    deltaY = dy;
    
    % Total number of unknowns
    N_total = (Nx - 1) * (Ny - 1);
    
    % Initialize sparse matrix
    A = sparse(N_total, N_total);
    
    % Loop to populate the matrix
    for i = 1:(Nx - 1)
        for j = 1:(Ny - 1)
            % Calculate global index for point (i, j)
            index = (j-1) * (Nx-1) + i;
            
            % Diagonal element
            n_value = n_matrix(i, j); % Refractive index at (i, j)
            A(index, index) = -2 / (deltaX^2) - 2 / (deltaY^2) + (k0^2) * (n_value^2);
            
            % Off-diagonal elements (vertical neighbors)
            if j > 1 % Not on the top boundary
                A(index, index - (Nx-1)) = 1 / (deltaX^2);
            end
            if j < Ny-1 % Not on the bottom boundary
                A(index, index + (Nx-1)) = 1 / (deltaX^2);
            end
            
            % Off-diagonal elements (horizontal neighbors)
            if i > 1 % Not on the left boundary
                A(index, index - 1) = 1 / (deltaY^2);
            end
            if i < Nx-1 % Not on the right boundary
                A(index, index + 1) = 1 / (deltaY^2);
            end
        end
    end
    
    % Display the resulting matrix
    % disp('Matrix A:');
    % disp(A);
    
    % Visualize the matrix structure
    % figure;
    % spy(A);
    % title('Structure of Matrix A');
    % xlabel('Column Index');
    % ylabel('Row Index');
    
    % Compute eigenvalues and eigenvectors
    disp('Computing eigenvalues and eigenvectors...');
    [eigenvectors, eigenvalues_matrix] = eigs(A, N_total); % Compute N_total eigenvalues/vectors
    eigenvalues = diag(eigenvalues_matrix); % Extract eigenvalues
    
    % Display eigenvalues
    % disp('Eigenvalues:');
    % disp(eigenvalues);
    
    % Ask the user to select an index
    % index = input(['Select an eigenvalue index (1 to ', num2str(N_total), '): ']);
    index = modeField.Value;

    % Ensure the index is valid
    while index < 1 || index > N_total
        disp('Invalid index. Please select a valid index.');
        index = input(['Select an eigenvalue index (1 to ', num2str(N_total), '): ']);
    end
    
    % Display the selected eigenvalue and corresponding eigenvector
    selected_eigenvalue = eigenvalues(index);
    selected_eigenvector = eigenvectors(:, index);
    
    % disp(['Selected Eigenvalue: ', num2str(selected_eigenvalue)]);
    % disp('Corresponding Eigenvector:');
    % disp(selected_eigenvector);
    
    % Plot the selected eigenvector as a mode shape
    figure;
    mode_shape = reshape(selected_eigenvector, [Nx-1, Ny-1]);
    imagesc(abs(mode_shape)');
    colorbar;
    title(['Eigenvector (Mode Shape) for Eigenvalue ', num2str(selected_eigenvalue)]);
    xlabel('x-axis');
    ylabel('y-axis');
    set(gca, 'YDir', 'normal'); % Ensure correct y-axis direction
    
    % -------------------- Confinement Factor Calculation -------------------- %
    % Compute the total energy (sum of |E|^2 over the entire domain)
    total_energy = sum(sum(abs(mode_shape).^2));
    
    % For a planar waveguide, define the waveguide region as:
    % x from 0 to overallLength, and y from substrateThickness to (substrateThickness + waveguideThickness)
    waveguide_mask = (x(1:end-1) >= 0 & x(1:end-1) <= overallLength.Value)' & ...
                     (y(1:end-1) >= substrateThickness.Value & y(1:end-1) <= (substrateThickness.Value + waveguideThickness.Value));
    
    % Compute the energy confined in the waveguide region
    waveguide_energy = sum(sum(abs(mode_shape(waveguide_mask)).^2));
    
    % Calculate the Confinement Factor
    confinement_factor = waveguide_energy / total_energy;
    
    % Display the Confinement Factor in the console
    disp(['Confinement Factor: ', num2str(confinement_factor)]);
    
    % Show the result in the final answer box in the GUI
    finalAnswerBox.Value = num2str(confinement_factor);

end

% Define the function to show inputs for Rectangular Waveguide
function showRectangularInputs()
    % Create a new figure for Rectangular inputs
    fig = uifigure('Position', [30, 50, 1200, 500], 'Name', 'Rectangular Waveguide Inputs');
    
    % Add an image (enlarged and positioned)
    axImage = uiaxes(fig, 'Position', [450, 90, 750, 400]); % Larger position
    imshow('./images/square_prosess.png', 'Parent', axImage, 'Border', 'tight'); % Remove border
    axImage.XColor = 'none';
    axImage.YColor = 'none';

    % Labels and input fields for refractive indices
    uilabel(fig, 'Position', [20, 450, 200, 20], 'Text', 'Refractive Index of Air:');
    airIndex = uieditfield(fig, 'numeric', 'Position', [280, 450, 100, 20], 'Value', 1);

    uilabel(fig, 'Position', [20, 420, 200, 20], 'Text', 'Refractive Index of Waveguide:');
    waveguideIndex = uieditfield(fig, 'numeric', 'Position', [280, 420, 100, 20], 'Value', 3.3);

    uilabel(fig, 'Position', [20, 390, 200, 20], 'Text', 'Refractive Index of Substrate:');
    substrateIndex = uieditfield(fig, 'numeric', 'Position', [280, 390, 100, 20], 'Value', 1.5);

    % Labels and input fields for overall dimensions
    uilabel(fig, 'Position', [20, 350, 250, 20], 'Text', 'Overall Dimensions (X Start, X End):');
    overallXStart = uieditfield(fig, 'numeric', 'Position', [280, 350, 50, 20], 'Value', 0);
    overallXEnd = uieditfield(fig, 'numeric', 'Position', [340, 350, 50, 20], 'Value', 10);

    uilabel(fig, 'Position', [20, 320, 250, 20], 'Text', 'Overall Dimensions (Y Start, Y End):');
    overallYStart = uieditfield(fig, 'numeric', 'Position', [280, 320, 50, 20], 'Value', 0);
    overallYEnd = uieditfield(fig, 'numeric', 'Position', [340, 320, 50, 20], 'Value', 10);

    % Labels and input fields for substrate dimensions
    uilabel(fig, 'Position', [20, 280, 250, 20], 'Text', 'Substrate Dimensions (X Start, X End):');
    substrateXStart = uieditfield(fig, 'numeric', 'Position', [280, 280, 50, 20], 'Value', 0);
    substrateXEnd = uieditfield(fig, 'numeric', 'Position', [340, 280, 50, 20], 'Value', 10);

    uilabel(fig, 'Position', [20, 250, 250, 20], 'Text', 'Substrate Dimensions (Y Start, Y End):');
    substrateYStart = uieditfield(fig, 'numeric', 'Position', [280, 250, 50, 20], 'Value', 0);
    substrateYEnd = uieditfield(fig, 'numeric', 'Position', [340, 250, 50, 20], 'Value', 6);

    % Labels and input fields for waveguide dimensions
    uilabel(fig, 'Position', [20, 210, 250, 20], 'Text', 'Waveguide Dimensions (X Start, X End):');
    waveguideXStart = uieditfield(fig, 'numeric', 'Position', [280, 210, 50, 20], 'Value', 2);
    waveguideXEnd = uieditfield(fig, 'numeric', 'Position', [340, 210, 50, 20], 'Value', 8);

    uilabel(fig, 'Position', [20, 180, 250, 20], 'Text', 'Waveguide Dimensions (Y Start, Y End):');
    waveguideYStart = uieditfield(fig, 'numeric', 'Position', [280, 180, 50, 20], 'Value', 6);
    waveguideYEnd = uieditfield(fig, 'numeric', 'Position', [340, 180, 50, 20], 'Value', 8);

    % Labels and input fields for mesh points
    uilabel(fig, 'Position', [20, 140, 200, 20], 'Text', 'Number of Mesh Points (X, Y):');
    meshPointsX = uieditfield(fig, 'numeric', 'Position', [280, 140, 50, 20], 'Value', 49);
    meshPointsY = uieditfield(fig, 'numeric', 'Position', [340, 140, 50, 20], 'Value', 49);

    % New fields for Mode and Lambda
    uilabel(fig, 'Position', [20, 100, 200, 20], 'Text', 'Mode:');
    modeField = uieditfield(fig, 'numeric', 'Position', [280, 100, 100, 20], 'Value', 1);

    uilabel(fig, 'Position', [20, 70, 200, 20], 'Text', 'Lambda (Wavelength):');
    lambdaField = uieditfield(fig, 'numeric', 'Position', [280, 70, 100, 20], 'Value', 1.55);

    % Create a styled label for the final answer box
    uilabel(fig, 'Position', [800, 100, 200, 30], 'Text', 'Confinement Factor', ...
        'FontSize', 14, 'FontWeight', 'bold', 'FontColor', [0 0 0.5]); % Dark blue color
    
    % Create a box for the final answer at the bottom right
    finalAnswerBox = uitextarea(fig, 'Position', [760, 70, 200, 30], 'Editable', 'on', ...
        'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', [0.9 0.9 0.9]); % Light gray background

    % Add a confirm button
    uibutton(fig, 'Position', [200, 20, 100, 30], 'Text', 'Submit', ...
    'ButtonPushedFcn', @(btn, event) submitRectangularInputs(airIndex, waveguideIndex, substrateIndex, ...
    overallXStart, overallXEnd, overallYStart, overallYEnd, ...
    substrateXStart, substrateXEnd, substrateYStart, substrateYEnd, ...
    waveguideXStart, waveguideXEnd, waveguideYStart, waveguideYEnd, ...
    meshPointsX, meshPointsY, modeField, lambdaField, finalAnswerBox));

end


function submitRectangularInputs(airIndex, waveguideIndex, substrateIndex, overallXStart, overallXEnd, overallYStart,...
    overallYEnd, substrateXStart, substrateXEnd, substrateYStart, substrateYEnd, waveguideXStart, waveguideXEnd, waveguideYStart,...
    waveguideYEnd, meshPointsX, meshPointsY, modeField, lambdaField, finalAnswerBox)

    % Gather input values
    inputs = struct( ...
        'AirIndex', airIndex.Value, ...
        'WaveguideIndex', waveguideIndex.Value, ...
        'SubstrateIndex', substrateIndex.Value, ...
        'OverallX', [overallXStart.Value, overallXEnd.Value], ...
        'OverallY', [overallYStart.Value, overallYEnd.Value], ...
        'SubstrateX', [substrateXStart.Value, substrateXEnd.Value], ...
        'SubstrateY', [substrateYStart.Value, substrateYEnd.Value], ...
        'WaveguideX', [waveguideXStart.Value, waveguideXEnd.Value], ...
        'WaveguideY', [waveguideYStart.Value, waveguideYEnd.Value], ...
        'MeshPoints', [meshPointsX.Value, meshPointsY.Value], ...
        'ModeField', modeField.Value, ...
        'LambdaField', lambdaField.Value ...
    );

    disp('Rectangular Wave-guide Inputs:');
    disp(inputs);
    % Proceed with further processing or computation

    % --- Your Computational Code Starts Here ---
    clc;
    close all;
    
    % Input values for domain and mesh size
    Lx = overallXEnd.Value - overallXStart.Value;
    Ly = overallYEnd.Value - overallYStart.Value;
    
    % Calculate the number of mesh points
    Nx = meshPointsX.Value;
    Ny = meshPointsY.Value;
    
    dx = Lx / Nx;
    dy = Ly / Ny;
    
    % Generate grid coordinates
    x = linspace(0, Lx, Nx + 1); % x-coordinates
    y = linspace(0, Ly, Ny + 1); % y-coordinates
    
    % Initialize the n(x, y) matrix with air (n = 1.0)
    n_matrix = ones(Nx - 1, Ny - 1);
    
    % Input substrate region
    substrate_x_min = substrateXStart.Value;
    substrate_x_max = substrateXEnd.Value;
    substrate_y_min = substrateYStart.Value;
    substrate_y_max = substrateYEnd.Value;
    
    % Assign refractive index for the substrate
    for i = 1:(Nx - 1)
        for j = 1:(Ny - 1)
            if x(i) >= substrate_x_min && x(i) <= substrate_x_max && ...
               y(j) >= substrate_y_min && y(j) <= substrate_y_max
                n_matrix(i, j) = substrateIndex.Value;
            end
        end
    end
    
    % Input waveguide parameters
    waveguide_type = 'Rectangular';
    
    if strcmpi(waveguide_type, 'Rectangular')
        % Base rectangle of the waveguide
        base_x_min = waveguideXStart.Value;
        base_x_max = waveguideXEnd.Value;
        base_y_min = waveguideYStart.Value;
        base_y_max = waveguideYEnd.Value;
    
        % Assign refractive index for the base rectangle
        for i = 1:(Nx - 1)
            for j = 1:(Ny - 1)
                if x(i) >= base_x_min && x(i) <= base_x_max && ...
                   y(j) >= base_y_min && y(j) <= base_y_max
                    n_matrix(i, j) = waveguideIndex.Value;
                end
            end
        end
    end
    
    % Visualize the refractive index distribution
    figure;
    imagesc(x, y, n_matrix');
    colorbar;
    title('n(x, y) Distribution');
    xlabel('x-axis');
    ylabel('y-axis');
    set(gca, 'YDir', 'normal');

    % Define wavelength and compute k0
    lambda = lambdaField.Value;
    k0 = 2 * pi / lambda;
    
    % Mesh sizes
    deltaX = dx;
    deltaY = dy;
    
    % Total number of unknowns
    N_total = (Nx - 1) * (Ny - 1);
    
    % Initialize sparse matrix A
    A = sparse(N_total, N_total);
    
    % Populate matrix A using finite difference approximation
    for i = 1:(Nx - 1)
        for j = 1:(Ny - 1)
            index = (j - 1) * (Nx - 1) + i;
            n_value = n_matrix(i, j);
            A(index, index) = -2 / (deltaX^2) - 2 / (deltaY^2) + (k0^2) * (n_value^2);
            if j > 1
                A(index, index - (Nx - 1)) = 1 / (deltaX^2);
            end
            if j < Ny - 1
                A(index, index + (Nx - 1)) = 1 / (deltaX^2);
            end
            if i > 1
                A(index, index - 1) = 1 / (deltaY^2);
            end
            if i < Nx - 1
                A(index, index + 1) = 1 / (deltaY^2);
            end
        end
    end
    
    % Compute eigenvalues and eigenvectors
    disp('Computing eigenvalues and eigenvectors...');
    [eigenvectors, eigenvalues_matrix] = eigs(A, N_total);
    eigenvalues = diag(eigenvalues_matrix);
    
    % Use modeField.Value as the index for selecting the eigenvalue/eigenvector
    index = modeField.Value;
    while index < 1 || index > N_total
        disp('Invalid index. Please select a valid index.');
        index = input(['Select an eigenvalue index (1 to ', num2str(N_total), '): ']);
    end
    
    selected_eigenvalue = eigenvalues(index);
    selected_eigenvector = eigenvectors(:, index);
    
    % Plot the selected eigenvector as a mode shape
    figure;
    mode_shape = reshape(selected_eigenvector, [Nx - 1, Ny - 1]);
    imagesc(abs(mode_shape)');
    colorbar;
    title(['Eigenvector (Mode Shape) for Eigenvalue ', num2str(selected_eigenvalue)]);
    xlabel('x-axis');
    ylabel('y-axis');
    set(gca, 'YDir', 'normal'); % Ensure correct y-axis direction

    % -------------------- Confinement Factor Calculation -------------------- %
    % Compute the total energy over the entire domain (sum of |E|^2)
    total_energy = sum(sum(abs(mode_shape).^2));
    
    % Create a mask for the waveguide region using the provided waveguide dimensions.
    % Note: Since mode_shape is defined on the internal grid, use x(1:end-1) and y(1:end-1)
    waveguide_mask = (x(1:end-1) >= waveguideXStart.Value & x(1:end-1) <= waveguideXEnd.Value)' & ...
                     (y(1:end-1) >= waveguideYStart.Value & y(1:end-1) <= waveguideYEnd.Value);
    
    % Compute the energy confined in the waveguide region
    waveguide_energy = sum(sum(abs(mode_shape(waveguide_mask)).^2));
    
    % Calculate the Confinement Factor
    confinement_factor = waveguide_energy / total_energy;
    
    % Display the Confinement Factor in the command window
    disp(['Confinement Factor: ', num2str(confinement_factor)]);
    
    % Show the result in the final answer box in the GUI
    finalAnswerBox.Value = num2str(confinement_factor);
end
