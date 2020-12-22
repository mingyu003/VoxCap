function [res_matvect]=lse_matvect_mult_charge_pecdiel_slow(ch_coef_dum,Z_mat)
% if fl_precond==2
%     for ll=1:size(ids,1)
%         ch_coef_dum(ids{ll},:)=precond{ll}*ch_coef_dum(ids{ll},:);
%     end
%     if isempty(ord)==0
%         for kk=1:size(ord,1)
%             ch_coef_dum(ord{kk},:)=(1/ss(end))*ch_coef_dum(ord{kk},:);
%         end
%     end
% end
res_matvect=Z_mat*ch_coef_dum;
fprintf ('.') ;
end