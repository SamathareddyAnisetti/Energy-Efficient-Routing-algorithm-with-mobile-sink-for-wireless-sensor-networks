clc
clear
load pegasis_data
load leach_data

%% Pegasis plot
    % Plotting Simulation Results "Operating Nodes per Transmission" %
    figure
    plot(1:rnd,op(1:rnd),'-r','Linewidth',2);
    title ({'PEGASIS'; 'Network Life Time';})
    xlabel 'Transmissions';
    ylabel 'Operational Nodes';
    hold on;
    
    
    % Plotting Simulation Results "Average Energy consumed by a Node per Transmission" %
    figure
    plot(1:rnd,abs(avg_energy),'-r','Linewidth',2);
    title ({'PEGASIS'; 'Average Energy consumed by a Node per Transmission';})
    xlabel 'Transmissions';
    ylabel 'Energy ( J )';
    hold on;
    %% Leach plot
        % Plotting Simulation Results "Operating Nodes per Transmission" %
    figure
    plot(1:r,AliveSensors(1:r),'-r','Linewidth',2);
    title ({'LEACH'; 'Network Life Time';})
    xlabel('Transmissions');
    ylabel('Operational Nodes');
    hold on;
    
    
    % Plotting Simulation Results "Average Energy consumed by a Node per Transmission" %
    figure
    plot(1:r,AvgEnergyAllSensor(2:r+1),'-r','Linewidth',2);
    title ({'LEACH'; 'Average Energy consumed by a Node per Transmission';})
    xlabel( 'Transmissions');
    ylabel( 'Energy ( J )');
    hold on;
    
    %% Comparison
    % Plotting Simulation Results "Operating Nodes per Transmission" %
    figure
    plot(1:rnd,op(1:rnd),'-r','Linewidth',2);
    hold on
    plot(1:r,AliveSensors(1:r),'-g','Linewidth',2);
    legend('PEGASIS','LEACH')
    title ({'Network Life Time';})
    xlabel 'Transmissions';
    ylabel 'Operational Nodes';
    

    
    % Plotting Simulation Results "Average Energy consumed by a Node per Transmission" %
    figure
    plot(1:rnd,abs(avg_energy),'-r','Linewidth',2);
    hold on
    plot(1:r,AvgEnergyAllSensor(2:r+1),'-g','Linewidth',2);
    legend('PEGASIS','LEACH')
    title ({'Average Energy consumed by a Node per Transmission';})
    xlabel 'Transmissions';
    ylabel 'Energy ( J )';
    hold on;