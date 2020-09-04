function [topLeftCornerTemp, bottomRightCornerTemp] = HeatTransferProject(m)
% creating constants
Tinf = 25; % T infinity
Kx = 5; % thermal conductivity in x
Ky = 10; % thermal conductivity in y
x = 25 * .001; % converting width into meters
y = 100 * .001; % converting height into meters
dx = x / (m-1); % finding dx
dy = y / (m-1); % finding dy
h1 = 100; % horizontal convection coefficient
qdot = 150000; % heat generation
coeff = []; % initializes an empty vector for the temperatures
A = []; % initializes an empty vector to be used for final matrix
h2 = [10 30 50 70 100]; % creating a vector for the vertical convection coefficient
i = 1:length(h2); % initializing the index to use with the different h2 values
Cvec = []; % initializes an empty vector for the non-solution solutions
j = 1; % initializes j at 1 
horzTemp1 = []; % initializes vector for first horizontal temp
horzTemp2 = []; % initializes vector for second horizontal temp
horzTemp3 = []; % initializes vector for third horizontal temp
horzTemp4 = []; % initializes vector for fourth horizontal temp
horzTemp5 = []; % initializes vector for fifth horizontal temp
vertTemp1 = []; % initializes vector for first vertical temp
vertTemp2 = []; % initializes vector for second vertical temp
vertTemp3 = []; % initializes vector for third vertical temp
vertTemp4 = []; % initializes vector for fourth vertical temp
vertTemp5 = []; % initializes vector for fifth vertical temp
topLeftCornerVec = []; % initializes vector for top left corner vec
bottomRightCornerVec = []; % initializes vector for bottom right corner vec
T_all = [];
for h2 = h2(i) % creating a for loop to loop through different h values
    j = 1; % re initializing j
    A = []; % re initializing A
    Cvec = []; % re initializing A
    while j <= m^2 % creates while loop for the length of j
        coeff = zeros(1,m^2); % creating a vector of length m^2
        if j == 1 % condition for top left point
            coeff(j) = (-Ky*(dx/dy)*(1/2))+(-Kx*(dy/dx)*(1/2))+(-h1*(dx/2)); % coefficient for top left point
            coeff(j+m) = (Ky*(dx/dy)*(1/2)); % coefficient for point below
            coeff(j+1) = (Kx*(dy/dx)*(1/2)); % coefficient for point to the right of current point
            Cvec = [Cvec; (-qdot*dx*dy*(1/4)+(-h1*Tinf*(dx/2)))]; % adds the current non-solution to the matrix of non-solutions
            j = j + 1; % increments j
            A = [A; coeff]; % adding coeff to the A matrix
        elseif j > 1 && j < m % condition for points between the top corners
            coeff(j) = (-Kx*(dy/dx)+(-Ky*(dx/dy))+(-h1*dx)); % coefficient for the point
            coeff(j-1) = (Kx*(dy/dx)*(1/2)); % coefficient for the point to the left 
            coeff(j+1) = (Kx*(dy/dx)*(1/2)); % coefficient for the point to the right
            coeff(j+m) = (Ky*(dx/dy)); % coefficient for the point below
            Cvec = [Cvec; (-qdot*((dx*dy)/2)-(h1*dx*Tinf))]; % adds another row to the Cvec
            j = j + 1; % increments j
            A = [A; coeff]; % adds coeff to the A matrix
        elseif j == m % condition for top right point
            coeff(j) = -Kx*(dy/dx)*(1/2)-(h2*(dy*(1/2)))-(Ky*(dx/dy)*(1/2))-(h1*(dx*(1/2)));
            coeff(j-1) = (Kx*(dy/dx)*(1/2)); % coefficient for point to the left
            coeff(j+m) = (Ky*(dx/dy)*(1/2)); % coefficient for the point below
            Cvec = [Cvec; (-qdot*((dx*dy)/4)-(Tinf*h1*dx*(1/2))-(Tinf*h2*(dy/2)))]; % adds a row to the Cvec matrix
            j = j + 1; % increments j
            A = [A; coeff]; % adds coeff to the A matrix
        elseif mod(j-1,m) == 0 && j ~= (m^2)-m+1 % condition for points on the left edge
            coeff(j) = ((-Ky*(dx/dy))-(Kx*(dy/dx))); % coeff for current point
            coeff(j+m) = (Ky*(dx/dy)*(1/2)); % coeff for point below
            coeff(j+1) = (Kx*(dy/dx)); % coeff for point to the right of the current point
            coeff(j-m) = (Ky*(dx/dy)*(1/2)); % coeff for point above
            j = j + 1; % increments j
            Cvec = [Cvec; (-qdot*(dx*dy)*(1/2))]; % adding a row to Cvec
            A = [A; coeff]; % adds coeff to the A matrix
        elseif j == m^2-m+1 % condition for bottom left point
            coeff(j) = (-Kx*(dy/dx)*(1/2))-(Ky*(dx/dy)*(1/2)); % coeff for current point
            coeff(j+1) = (Kx*(dy/dx)*(1/2)); % coeff for point to the left
            coeff(j-m) = (Ky*(dx/dy)*(1/2)); % coeff for point above
            A = [A; coeff]; % add row to A matrix
            Cvec = [Cvec; ((-qdot*dx*dy)/4)]; % add row to Cvec
            j = j + 1; % increments j
        elseif j == m^2 % bottom right corner
            coeff(j) = (-Kx*(dy/dx)*(1/2))-(h2*dy*(1/2))-(Ky*(dx/dy)*(1/2)); % coeff for current point
            coeff(j-1) = (Kx*(dy/dx)*(1/2)); % coeff for point to the the left
            coeff(j-m) = (Ky*(dx/dy)*(1/2)); % coeff for point above
            Cvec = [Cvec; ((-qdot*dx*dy*(1/4))-(h2*Tinf*dy*(1/2)))]; %add a row to Cvec
            A = [A; coeff]; % add row to A matrix
            j = j + 1; % increment j
        elseif mod(j,m) == 0 % right edge
            coeff(j) = (-Kx*(dy/dx))-(Ky*(dx/dy))-h2*dy; % coeff for current point
            coeff(j-1) = (Kx*(dy/dx)); % coeff for point to the left
            coeff(j+m) = (Ky*(dx/dy)*(1/2)); % coeff for point below
            coeff(j-m) = (Ky*(dx/dy)*(1/2)); % coeff for point above
            A = [A; coeff]; % add row to A matrix
            Cvec = [Cvec; (-qdot*(dx*dy)*(1/2))-(h2*(dy*Tinf))]; % add a row to Cvec
            j = j + 1; % increments j
        elseif j > m^2-m+1 && j < m^2 % condition for bottom edge
            coeff(j) = (-Kx*(dy/dx))-(Ky*(dx/dy)); % coeff for current point
            coeff(j-1) = (Kx*(dy/dx)*(1/2)); % coeff for point to the left
            coeff(j+1) = (Kx*(dy/dx)*(1/2)); % coeff for poiont to the right
            coeff(j-m) = (Ky*(dx/dy)); % coeff for point above
            Cvec = [Cvec; (-qdot*dx*dy*(1/2))]; % add row to Cvec
            A = [A; coeff]; % add row to A matrix
            j = j + 1; % increment 1
        else
            coeff(j) = (-2*(Kx*dy)/dx)-(2*(Ky*dx)/dy); % coeff at current point
            coeff(j-1) = (Kx*(dy/dx)); % coeff for point to the left
            coeff(j+1) = (Kx*(dy/dx)); % coeff for point to the right
            coeff(j+m) = (Ky*(dx/dy)); % coeff for point below
            coeff(j-m) = (Ky*(dx/dy)); % coeff for point above
            Cvec = [Cvec; (-qdot*dx*dy)]; % add row to Cvec
            A = [A; coeff]; % add row to the A matrix
            j = j + 1; % increments j
        end
    end
    % solve for temperatures
    T = (A^-1)*Cvec; % solution matrix
    a=((m-1)/2)+1; % vertical temperatures
    b=((m-1)/4)*m; % how many temperatures to find
    c=(m^2)-((m-1)/2); % last temperature on bottom row
    ind = linspace(a,c,5); % 5 evenly spaced points
    T_vert=T(ind); % Vertical temperatures along center line
    d=b/m; % used to find the temperature
    horzind=ind(3); % middle vertical value which lines up with horizontal temperatures
    horzind = linspace(horzind-2*d,horzind+2*d,5); % evenly spaced on horizontal
    T_horz=T(horzind); % horizontal temperatures along center line
    horzTemp1 = [horzTemp1 T_horz(1)]; % fills in the vector
    horzTemp2 = [horzTemp2 T_horz(2)]; % fills in the vector
    horzTemp3 = [horzTemp3 T_horz(3)]; % fills in the vector
    horzTemp4 = [horzTemp4 T_horz(4)]; % fills in the vector
    horzTemp5 = [horzTemp5 T_horz(5)]; % fills in the vector
    vertTemp1 = [vertTemp1 T_vert(1)]; % fills in the vector
    vertTemp2 = [vertTemp2 T_vert(2)]; % fills in the vector
    vertTemp3 = [vertTemp3 T_vert(3)]; % fills in the vector
    vertTemp4 = [vertTemp4 T_vert(4)]; % fills in the vector
    vertTemp5 = [vertTemp5 T_vert(5)]; % fills in the vector
    topLeftCornerVec = [topLeftCornerVec T(1)];
    bottomRightCornerVec = [bottomRightCornerVec T(m^2)];
    T_all = [T_all T];
end
    topLeftCornerTemp = topLeftCornerVec(1); % finds top left corner temp for last h2
    bottomRightCornerTemp = bottomRightCornerVec(1); % finds bottom right corner temp for last h2
    h2 = [10 30 50 70 100]; % restating the h coefficients
    figure(1); % initializes the first figure
    plot(h2,horzTemp1,'-bo',h2,horzTemp2,'-r*',h2,horzTemp3,'-k+',h2,horzTemp4,'-gd',h2,horzTemp5,'-ms'); % plots
    xlabel('h (W/m^2k)'); % creates label for x axis
    ylabel('Temperature (C)'); % creates label for y axis
    xlim([0 110]); % setting x axis range
    title('Horizontal Temperature Solutions'); % creates a title
    legend('1','2','3','4','5'); % Creates a legend
    figure(2); % initializes the second figure
    plot(h2,vertTemp5,'-bo',h2,vertTemp4,'-r*',h2,vertTemp3,'-k+',h2,vertTemp2,'-gd',h2,vertTemp1,'-ms'); % plots)
    xlabel('h (W/m^2k)'); % creates label for x axis
    ylabel('Temperature (C)'); % creates label for y axis
    xlim([0 110])
    title('Vertical Temperature Solutions'); % creates a title
    legend('1','2','3','4','5'); % Creates a legend
    
    T_h10 = T_all(:,1);
    T_contour_h10 = flip(transpose(reshape(T_h10,[m,m])));
    subplot(1,2,1)
    contourf(T_contour_h10,100,'LineColor',"none")
    colormap("autumn")
    colorbar
    xticks(linspace(0,m,6))
    yticks(linspace(0,m,6))
    xticklabels({'0','5','10','15','20','25'})
    yticklabels({'0','20','40','60','80','100'})
    set(gcf,'position',[100,100,350,500])
    xlabel('width [mm]')
    ylabel('height [mm]')
    title('Convection on right side: 10 [W/m^2-K]')
    print -r500
    
    T_h100 = T_all(:,5);
    T_contour_h100 = flip(transpose(reshape(T_h100,[m,m])));
    subplot(1,2,2)
    contourf(T_contour_h100,100,'LineColor',"none")
    colormap("autumn")
    colorbar
    xticks(linspace(0,m,6))
    yticks(linspace(0,m,6))
    xticklabels({'0','5','10','15','20','25'})
    yticklabels({'0','20','40','60','80','100'})
    set(gcf,'position',[100,100,350,500])
    xlabel('width [mm]')
    ylabel('height [mm]')
    title('Convection on right side: 100 [W/m^2-K]')
    print -r500
    
end
