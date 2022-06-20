function basicwaitbar(i,n,msg)
%fwait = waitbar(0,msg);

waitbar(i/n,strcat(msg," ",sprintf('%g%%',round(i/n*100,2))))

%close(f)
end