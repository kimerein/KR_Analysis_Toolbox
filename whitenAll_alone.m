function output=whitenAll_alone(input,isCx,isHigh)

if isCx==false
    for i=1:size(input,3)
        temp=reshape(input(:,:,i),size(input,1),size(input,2));
        temp(:,1:3)=temp(:,1:3)./(1.9*0.17.*(0.7225/0.1121));
        temp(:,4)=temp(:,4)./(1.4*0.73.*(0.7225/0.5035));
        temp(:,2)=1.07*temp(:,2);
        temp(:,4)=1.07*temp(:,4);
        temp(:,3)=1.0475*temp(:,3);
        temp(:,1:3)=temp(:,1:3)+0.1;
        temp(:,4)=temp(:,4)+0.1;
        leaveMinAt=nanmin(temp(:,3));
        temp(:,3)=temp(:,3)*1.5-1.5*leaveMinAt-0.1;
        temp(:,2)=temp(:,2)*1.2-0.9*leaveMinAt+0.18;
        temp(:,4)=temp(:,4)*1.2-0.9*leaveMinAt+0.18;
        output(:,:,i)=temp;
    end
else
    for i=1:size(input,3)
%         temp=reshape(input(:,:,i),size(input,1),size(input,2));
%         temp(:,1:3)=temp(:,1:3)./((3.3*0.17.*(0.7225/0.1121))*0.37);
%         temp(:,4)=temp(:,4)./((2.0*0.73.*(0.7225/0.5035))*0.39);
%         temp(:,2)=1.07*temp(:,2);
%         temp(:,1)=1.05*temp(:,1);
%         temp(:,3)=1.06*temp(:,3);
%         temp(:,4)=0.75*temp(:,4);
%         temp=temp-repmat(nanmin(temp,[],2),1,size(temp,2));
%         temp=temp./repmat(nanmax(temp,[],2),1,size(temp,2));
%         output(:,:,i)=temp;
         if isHigh==0
            temp=reshape(input(:,:,i),size(input,1),size(input,2));
            temp(:,1:3)=temp(:,1:3)./((3.3*0.17.*(0.7225/0.1121))*1*1.1);
            temp(:,4)=temp(:,4)./((2.0*0.73.*(0.7225/0.5035))*0.7*1.1);
            temp(:,2)=0.95*temp(:,2);
            temp(:,1)=1.0*temp(:,1);
            temp(:,3)=1.07*temp(:,3);
            temp(:,4)=0.42*temp(:,4);
            temp(:,5)=0.48*temp(:,5);
            temp(:,6)=0.95*temp(:,6);
            output(:,:,i)=temp;
            
         else
%             temp=reshape(input(:,:,i),size(input,1),size(input,2));
%             temp(:,1:3)=temp(:,1:3)./((3.3*0.17.*(0.7225/0.1121))*1.5);
%             temp(:,4)=temp(:,4)./((2.0*0.73.*(0.7225/0.5035))*0.9);
%             temp(:,2)=1.07*temp(:,2);
%             temp(:,1)=1.03*temp(:,1);
%             temp(:,3)=1.06*temp(:,3);
%             temp(:,4)=0.75*temp(:,4);
%             output(:,:,i)=temp;
            
            temp=reshape(input(:,:,i),size(input,1),size(input,2));
            temp(:,1:3)=temp(:,1:3)./((3.3*0.17.*(0.7225/0.1121))*1*1.3);
            temp(:,4)=temp(:,4)./((2.0*0.73.*(0.7225/0.5035))*0.7*1.3);
            temp(:,2)=0.95*temp(:,2);
            temp(:,1)=1.0*temp(:,1);
            temp(:,3)=1.07*temp(:,3);
            temp(:,4)=0.42*temp(:,4);
            temp(:,5)=0.43*temp(:,5);
            temp(:,6)=0.95*temp(:,6);
            output(:,:,i)=temp;
         end
    end
end

end