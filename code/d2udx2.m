function ut = d2udx2(t,u)
	%{
	Universidade Federal da Paraíba
	Centro de Tecnologia
	Departamento de Engenharia Mecânica

	Trabalho de Conclusão de Curso
	
	Análise Numérica de Condução Transiente com Termo Fonte Variável em
	Varetas Combustíveis de Reatores Nucleares pelo Método das Linhas 

	Discente: Thiago Ney Evaristo Rodrigues
	Orientador: Dr. Jacques César dos Santos
	%}

	% Problem parameters

    stcase = 1; % Case in study
    
	r0 = 0.001;
	rl = 1;
	nr = length(u);
    Bi = 15;
    Gast = 32.4;
    r = linspace(r0,rl,nr);
    
	% Partial derivative

	dr = (rl - r0)/(nr-1);
    dr2 = dr^2;
    
    ur = zeros(nr,1); % Preallocation
	urr = zeros(nr,1); % Preallocation

    % Partial derivative
    
	for i = 2:(nr-1)
		
        ur(i) = (u(i) - u(i-1))/dr;
        urr(i) = (u(i+1) - 2*u(i) + u(i-1))/dr2;
		
	end

	% Boundary condition

    ur(1) = (u(1) - u(2))/dr;
	urr(1) = 2*(u(2) - u(1))/dr2;
    ur(end) = (u(end) - u(end - 1))/dr;
	urr(end) = 2*(u(end - 1) - (Bi*dr + 1)*u(end))/dr2;

    % Partial differential equation

    ut = zeros(nr,1);
    
    if stcase == 1
            
        Glin = Gast*t;
        
        for i = 1:nr
            
            ut(i) = (1/r(i))*ur(i) + urr(i) + Glin;
        
        end

    elseif stcase == 2

        c3 = 4;
        Glin = Gast*(r.^2)*exp(c3*t);
        
        for i = 1:nr
            
            ut(i) = (1/r(i))*ur(i) + urr(i) + Glin(i);
            
        end

    elseif stcase == 3

        c1 = 1;
        Glin = Gast*(1 + c1*t);
        
        for i = 1:nr
            
            ut(i) = (1/r(i))*ur(i) + urr(i) + Glin;
        
        end

    else
        
        c2 = 1;
        c3 = 1;
        Glin = Gast*(1 + c2*(r.^2))*exp(c3*t);
        
        for i = 1:nr
            
            ut(i) = (1/r(i))*ur(i) + urr(i) + Glin(i);
            
        end

    end
  
end