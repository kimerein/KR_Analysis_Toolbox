function isRunning=plotRunning_wrapper(encoder,LED)

isRunning=nan(1,size(encoder,1));
for i=1:size(encoder,1)
    isRunning(i)=plotRunning(encoder,LED,i);
end