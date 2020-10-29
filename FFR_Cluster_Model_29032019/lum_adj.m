function [JtNa,JtK,JCl,Qtot,Nal,Kl,Cll,QwNa,QwK,QwCl] = ...
                              lum_adj(Nal,Kl,Cll,JtNa,JtK,JCl,Qtot,par)
 % #########################################################################
% -------------------------------------------------------------------------
% Author: Elias Siguenza
% Location: The University of Auckland, New Zealand
% Date: 22 March 2019
% Version: 1.1
% -------------------------------------------------------------------------
% Purpose:
% This function uses the adjacency matrix to calculate the luminal structure
% equations. 
% This is a dirty function and should be re-written.
% -------------------------------------------------------------------------
% #########################################################################
% Add upper to lower triangles:
    Qtot = (Qtot - triu(Qtot)) + triu(Qtot)';
    JtNa  = (JtNa - triu(JtNa)) + triu(JtNa)';
    JtK   = (JtK - triu(JtK)  ) + triu(JtK)';
    JCl  = (JCl - triu(JCl)  ) + triu(JCl)';                    
                          
    Nal = triu(Nal)';
    Kl = triu(Kl)';
    Cll = triu(Cll)';  
    
    [qa,qb]=find(Qtot');
    qa = [qb,qa];
    
    
    JtNad = zeros(19,1);
    JtKd = zeros(19,1);
    JCld = zeros(19,1);
    Qtotd = zeros(19,1);
    Nald = zeros(19,1);
    Kld = zeros(19,1);
    Clld = zeros(19,1);
    
    for j = 1:19
        Nald(j,1) = Nal(qa(j,1),qa(j,2));
        Kld(j,1)  = Kl(qa(j,1),qa(j,2));
        Clld(j,1) = Cll(qa(j,1),qa(j,2));
        Qtotd(j,1) = Qtot(qa(j,1),qa(j,2));
        JtNad(j,1) = JtNa(qa(j,1),qa(j,2));
        JtKd(j,1)  = JtK(qa(j,1),qa(j,2));
        JCld(j,1)  = JCl(qa(j,1),qa(j,2));
    end
    
    % Multiply by adjacency matrix
    JtNa = par.Adjs.*JtNad;
    JtK = par.Adjs.*JtKd;
    JCl = par.Adjs.*JCld;
    Qtot = par.Adjs.*Qtotd;
    Nal = par.Adjs.*Nald;
    Kl = par.Adjs.*Kld;
    Cll = par.Adjs.*Clld;
    
    
    % Precalculate the water/ion influx:
    QwNa = Qtot.*Nal;
    QwNa = QwNa - diag(diag(QwNa));
    QwK = Qtot.*Kl;
    QwK = QwK - diag(diag(QwK));
    QwCl = Qtot.*Cll;
    QwCl = QwCl - diag(diag(QwCl));
                          
    

    
    
end
    