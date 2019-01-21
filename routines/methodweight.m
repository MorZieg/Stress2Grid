%% method index, only used when you want to discriminate datas from different method
function ind_method = methodweight(methodname,fms,fmf,bo,dif,hf,gf,gv,oc,rest)
[m_,n_]=size(methodname);
ind_method = zeros(m_,n_);
for i_ = 1:m_
    for j_ = 1:n_
        methodij = char(methodname(i_,j_));
        methodij = methodij(1:2);
        switch methodij
            case 'OC'
                ind_method(i_,j_) = oc;
            case 'HF'
                ind_method(i_,j_) = hf;
            case 'BO'
                ind_method(i_,j_) = bo;
            case 'DI'
                ind_method(i_,j_) = dif;
            case 'FMS'
                ind_method(i_,j_) = fms;
            case 'FMF'
                ind_method(i_,j_) = fmf;
            case 'GF'
                ind_method(i_,j_) = gf;
            case 'GV'
                ind_method(i_,j_) = gv;
            otherwise
                ind_method(i_,j_) = rest;
        end
    end
end
end