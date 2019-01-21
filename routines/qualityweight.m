%% quality weight function
function quality_weight = qualityweight(data_quality,qwA,qwB,qwC,qwD,qwE)
[m_,n_]=size(data_quality);
quality_weight = zeros(m_,n_);
for i_ = 1:m_
    for j_ = 1:n_
        dqij = char(data_quality(i_,j_));
        switch dqij
            case 'A'
                quality_weight(i_,j_) = qwA;
            case 'B'
                quality_weight(i_,j_) = qwB;
            case 'C'
                quality_weight(i_,j_) = qwC;
            case 'D'
                quality_weight(i_,j_) = qwD;
            otherwise
                quality_weight(i_,j_) = qwE;
        end
    end
end
end