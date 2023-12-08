function [Error_b_split1,Error_m_split1] =calspliterror(train_result_matrix,numbenign,nummalignant)
b_col=size((train_result_matrix),2)-1;
m_col=size(train_result_matrix,2);
%%%%%%% for split 1 as b results in 4th col %%%%%%%%%%%%%
train_result_matrix;
Error_b_split1=0;
Error_m_split1=0;

for i= 1:size(train_result_matrix,1)

    if train_result_matrix(i,1)==2 && train_result_matrix(i,b_col)~=2
    Error_b_split1= Error_b_split1+1;
    end
    if train_result_matrix(i,1)==4 && train_result_matrix(i,b_col)~=0
    Error_b_split1= Error_b_split1+1;
    end

    if train_result_matrix(i,1)==2 && train_result_matrix(i,m_col)~=0
    Error_m_split1= Error_m_split1+1;
    end
    if train_result_matrix(i,1)==4 && train_result_matrix(i,m_col)~=4
    Error_m_split1= Error_m_split1+1;
    end
end
    num_splitb=sum(train_result_matrix(:,1)==2);
    num_splitm=sum(train_result_matrix(:,1)==4);
    numbenign;
    nummalignant;
    Error_b_split1=Error_b_split1/num_splitb;  %%%MAE
    Error_m_split1=Error_m_split1/num_splitm; 

end
