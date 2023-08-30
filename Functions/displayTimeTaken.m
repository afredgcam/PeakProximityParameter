function displayTimeTaken(timeTaken)

if timeTaken >= 60
    minutes = floor(timeTaken/60);
    secs = timeTaken - 60*minutes;
    disp(['time taken: ',num2str(minutes),' minutes and ',num2str(secs,2),' seconds'])
else
    disp(['time taken: ',num2str(timeTaken,3),' seconds'])
end

end