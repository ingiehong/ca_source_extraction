

% Plot spatial component and S,C
unit=reshape(P.sn,d1,d2)/max(P.sn); %zeros(d1,d2);
max_unit=reshape(P.sn,d1,d2)/max(P.sn); %zeros(d1,d2);

% Open decent size figure
figure
scrsz = get(groot,'ScreenSize');
set(gcf, 'Position',[50 50 scrsz(3)-100 scrsz(4)*2/3-150])
pan2 = panel();
pan2.pack('h', {1/4 1/4 []} )

for i=1:size(S_df,1)
    % Spatial component
    pan2(1).select()
    unit=reshape(A_or(:,i), d1, d2);
    imagesc(unit);
    axis image; axis equal; axis tight; axis ij;  
    title(['Unit ' num2str(i) ' / ' num2str(size(S_df,1))])
    
    % Maximal activation of Spatial component
    pan2(2).select()
    [~,max_i]=max(C_df(i,:));
    max_unit=max_unit + reshape( A_or(:,i).* Yr(:,max_i)/max(A_or(:,i).* Yr(:,max_i)), d1, d2);
    imagesc(max_unit);
    axis image; axis equal; axis tight; axis ij; 
    title(['Unit ' num2str(i) ' / ' num2str(size(S_df,1))])
    
    % S and C
    pan2(3).select()
    plot([S_df(i,:); C_df(i,:)]')
    pause
end