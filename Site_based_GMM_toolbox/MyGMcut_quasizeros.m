function [eq_cut,cutID,L_cut] = MyGMcut_quasizeros(cut_p_lower,cut_p_upper,dteq,eq)

if isempty(dteq)
    eq_cut=[];
    cutID=[];
    L_cut=[];
else
    L = numel(eq);
    teq = (0:L-1)*dteq;

    AI = cumtrapz(teq,eq.^2);
    AI_n = AI/AI(end)*100;
    factor = 100/(cut_p_upper-cut_p_lower);
    id1 = find(AI_n>=cut_p_lower,1,'first');
    id2 = find(AI_n>=cut_p_upper,1,'first');
    cutID = [id1,id2];
    eq_cut = eq(id1:id2)*sqrt(factor);
    L_cut  = numel(eq_cut);
end

end

%%

