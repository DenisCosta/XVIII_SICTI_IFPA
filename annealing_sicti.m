clear; clc; close all;
fprintf('\n===============================================================\n')
fprintf('\n  Instituto Federal de Educação, Ciência e Tecnologia do Pará')
fprintf('\n                  IFPA Câmpus Ananindeua')
fprintf('\n   Grupo de pesquisa: Gradiente de Modelagem Matemática e')
fprintf('\n               Simulação Computacional - GM²SC')
fprintf('\n    Código-fonte utilizado nas simulações da pesquisa:')
fprintf('\n       "Simulated Annealing aplicado na Otimização ')
fprintf('\n              Global de Funções Multimodais"')
fprintf('\n          Artigo apresentado no XVIII SICTI e XI SIMIT')
fprintf('\n Autores: Dr. Denis C. L. Costa, Ana Elisa Oliveira Saavedra,')
fprintf('\n Adamilie P. de Paula, Vinicius A. Barbosa, Maria Alice da Rocha,\n')
fprintf('\n===============================================================\n\n')
% Para executar em OCTAVE carregar o package statistics:
% pkg load statistics
% ==========================================================
%        Configuração dos Interruptores de Restrições
% ==========================================================
% Interruptor das Restrições: 1 => Ativo; 0 => Inativo
% Para a Função de Ackley padrão, manteve-se as restrições inativas (0).
global usar_h1; usar_h1 = 0;
global usar_h2; usar_h2 = 0;
global usar_g1; usar_g1 = 0;
global usar_g2; usar_g2 = 0;

% ==========================================================
%        Definição da Função Objetivo e as Penalidades
% ==========================================================
function z = objective_function(x)
    % Parâmetros padrão da Função de Ackley
    a = 20; b = 0.2; c = 2*pi;
    d = length(x); % Número de variáveis: d = 2
    
    soma_sq = x(1)^2 + x(2)^2;
    soma_cos = cos(c*x(1)) + cos(c*x(2));
    
    term1 = -a * exp(-b * sqrt(soma_sq / d));
    term2 = -exp(soma_cos / d);
    
    z = term1 + term2 + a + exp(1);
end
%% --------------------------------------------------------------------
function p = penalty_function(x)
    global usar_h1 usar_h2 usar_g1 usar_g2;
    p = 0;
    peso = 1000;
    % Restrições desativadas para busca global pura na Ackley
    if (usar_h1 == 1), p = p + peso * abs((x(1)^2 + x(2)^2) - 1); end
    if (usar_h2 == 1), p = p + peso * abs((x(1) + x(2)) - 2); end


% --- Restrições de Desigualdade (g) ---
    if (usar_g1 == 1)
        % Exemplo: x1 + x2 >= 1  => Violada se (x1 + x2) < 1
        if (x(1) + x(2) < 1)
            p = p + peso * abs(1 - (x(1) + x(2)));
        end
    end

    if (usar_g2 == 1)
        % Exemplo: x1^2 + x2 <= 0.5 => Violada se (x1^2 + x2) > 0.5
        if (x(1)^2 + x(2) > 0.5)
            p = p + peso * abs((x(1)^2 + x(2)) - 0.5);
        end
    end
end
%% --------------------------------------------------------------------

function f = fitness(x)
    f = objective_function(x) + penalty_function(x);
end

% ==========================================================
%         Parâmetros e Loop do Simulated Annealing 
% ==========================================================
T0 = 100;           % Temperatura inicial alta para explorar picos
T_min = 1e-6;       % Temperatura final mais baixa para precisão
alpha = 0.95;       % Fator de resfriamento (ajustado para 0.95 das iterações)
max_iter = 100;     % Iterações por nível de temperatura

% Ponto inicial aleatório dentro do domínio típico [-32, 32]
x_current = -32 + (32 - (-32)) * rand(1, 2);
f_current = fitness(x_current);

x_best = x_current;
f_best = f_current;

history = [];
T = T0;

while T > T_min
    for i = 1:max_iter
     % Perturbação: O passo deve ser proporcional ao domínio
     % Passo diminui levemente com a temperatura para refinamento
        x_new = x_current + randn(1, 2) * (T/T0 * 5 + 0.1); 
        
        % Garantir que os limites do domínio não sejam ultrapassados
        x_new = max(min(x_new, 32.768), -32.768);
        
        f_new = fitness(x_new);
        delta = f_new - f_current;
        
        % Critério de Metrópolis
        if (delta < 0) || (rand() < exp(-delta / T))
            x_current = x_new;
            f_current = f_new;
            
            if f_current < f_best
                x_best = x_current;
                f_best = f_current;
            end
        end
    end
    history = [history, f_best];
    T = T * alpha;
end

% ==========================================================
%         Representação Gráfica e dos Resultados
% ==========================================================
figure;
subplot(1,1,1);
semilogy(history, 'LineWidth', 3, 'Color', [1 0 0]);
grid on;
title('Convergência na Função de Ackley (Escala Log)');
xlabel('Evolução da Temperatura');
ylabel('Melhor Fitness Encontrado');

fprintf('    --- Resultados da Otimização ---\n');
fprintf('Ponto Ótimo Encontrado (x,y): [%.6f, %.6f]\n', x_best(1), x_best(2));
fprintf('Valor da Função f(x,y): %.6e\n', f_best);
fprintf('\n  [>>>>> Fim da Simulação <<<<<]\n')