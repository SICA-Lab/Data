% convert 120us bcg to 24us bcg using linear interpolation
function meas_bcg_new = BcgProcess(meas_bcg)
meas_bcg_new = cell(size(meas_bcg));
for ii = 1:27
    for jj = 1:27
        data = meas_bcg{jj,ii}(711:951,:);
        data1 = zeros(5,240);
        dd1 = data(1:end-1,2);
        dd2 = data(2:end,2);
        data1(1,:) = dd1;
        data1(2,:) = 0.8*dd1 + 0.2*dd2;
        data1(3,:) = 0.6*dd1 + 0.4*dd2;
        data1(4,:) = 0.4*dd1 + 0.6*dd2;
        data1(5,:) = 0.2*dd1 + 0.8*dd2;
        T = [67.7:0.02:91.68]'*1e-6;
        meas_bcg_new{jj,ii} = [T,reshape(data1*10,[],1)];
    end
end