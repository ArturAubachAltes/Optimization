clear all;
clc;

% Configuración inicial
modelName = 'Kalman_filter_LQG';
load_system(modelName);
set_param(modelName, 'AlgebraicLoopSolver', 'TrustRegion');

% Parámetros del sistema (reemplaza A y B por tus matrices específicas)
r = 0.006;
M_c = 0.135;
I = 0.0007176;
l = 0.2;
g = 9.81;
b = 0.00007892;
L = 0.046;
Rm = 12.5;
kb = 0.031;
kt = 0.031;
c = 0.63;
m = 0.1;
M = 0.136;
Er = 2*m*g*l;
n = 3;
AA = I*(M+m) + M*m*(l^2);
aa = (((m*l)^2)*g)/AA;
bb = ((I +m*(l^2))/AA)*(c + (kb*kt)/(Rm*(r^2)));
cc = (b*m*l)/AA;
dd = (m*g*l*(M+m))/AA;
ee = ((m*l)/AA)*(c + (kb*kt)/(Rm*(r^2)));
ff = ((M+m)*b)/AA;
mm = ((I +m*(l^2))*kt)/(AA*Rm*r);
nn = (m*l*kt)/(AA*Rm*r);
A = [0 0 1 0; 0 0 0 1; 0 aa -bb -cc; 0 dd -ee -ff];
B = [0; 0; mm; nn];
C=eye(4);
D=[0;0;0;0];

est_values=eig(A);
disp(est_values)
%Nuestro sistema es completamente controlable.
O = obsv(A, C);
rango_obs=rank(O);
disp(rango_obs);
%Nuestro sistema es completamente observable.
% Verificar controlabilidad del sistema
ctrlMat = ctrb(A, B);
rango=rank(ctrlMat);
disp(rango);
if rango < size(A, 1)
    error('El sistema no es controlable.');
end

% Estado deseado
desiredState = [0, pi, 0, 0];  % Modifica según tus requerimientos

% Valores iniciales de Q y R para la optimización
initialParams = [300, 500, 0, 0, 0.0035]; % Inicial Q y R
lb = [0, 0, 0, 0, 0.00035]; % Límites inferiores
ub = [3000, 3000, 0, 0, 5]; % Límites superiores

% Configuración de la optimización
options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...               % Algoritmo de optimización
    'Display', 'iter', ...                % Mostrar progreso iterativo
    'MaxFunctionEvaluations', 1000, ...   % Máximas evaluaciones de función
    'FiniteDifferenceStepSize', 100, ...    % Tamaño de paso de diferencia finita aumentado
    'StepTolerance', 0.01, ...            % Tolerancia de paso reducida
    'OptimalityTolerance', 0.01);         % Tolerancia de optimalidad reducida

% Llamada a fmincon para optimizar Q y R
optimalParams = fmincon(@(params) lqrCostFunction(params, A, B, modelName, desiredState), initialParams, [], [], [], [], lb, ub, [], options);

% Valores óptimos de Q y R
Q_optimal = diag(optimalParams(1:4));
R_optimal = optimalParams(5);

% Mostrar los resultados óptimos
disp('Valores óptimos de Q:');
disp(Q_optimal);
disp('Valor óptimo de R:');
disp(R_optimal);

% Cerrar el modelo de Simulink al finalizar
close_system(modelName, 0);

% Definir la función de costo
function cost = lqrCostFunction(params, A, B, modelName, desiredState)
    try
        % Extraer Q y R de los parámetros
        Q = diag(params(1:4));
        R = params(5);

        % Verificar que Q y R sean adecuadas
        if any(diag(Q) < 0)
            error('Todos los elementos de Q deben ser no negativos.');
        end
        if R <= 0
            error('R debe ser positiva.');
        end

        % Calcular la ganancia LQR K
        K = lqr(A, B, Q, R);

        % Asignar K al espacio de trabajo de MATLAB para que Simulink lo use
        assignin('base', 'KK', K);

        % Configurar el tiempo de simulación
        simulationTime = 30; % Duración de la simulación en segundos

        % Ejecutar la simulación
        simOut = sim(modelName, 'SimulationMode', 'normal', 'SaveOutput', 'on', 'SaveTime', 'on', 'StopTime', num2str(simulationTime));

        % Verificar las señales disponibles en la salida de simulación
        signalNames = simOut.who;
        disp('Señales disponibles en la salida de simulación:');
        disp(signalNames);

        % Extraer los datos de las señales correctas
        % Actualiza estos nombres según las señales correctas en tu modelo
        realState = simOut.get('realState'); % Ajusta según tu modelo
        desiredStateSim = simOut.get('desiredState'); % Ajusta según tu modelo

        % Asegurarse de que los datos sean matrices
        realState = realState.Data;
        desiredStateSim = desiredStateSim.Data;

        % Extraer los valores al final de la simulación
        finalRealState = realState(end, :); % Último estado real
        finalDesiredState = desiredStateSim(end, :); % Último estado deseado

        % Calcular el costo como la diferencia entre el estado real y el deseado
        cost = norm(finalRealState - finalDesiredState); % Usa la norma para evaluar la diferencia
        
    catch ME
        % En caso de error, mostrar el mensaje y devolver un costo alto
        disp(['Error en lqrCostFunction: ', ME.message]);
        cost = 1e10; % Valor alto para indicar fallo en la evaluación
    end
end
