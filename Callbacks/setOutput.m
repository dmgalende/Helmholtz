function setOutput(msje,handle,option)

if isempty(handle)
    disp(msje)
    return
elseif nargin == 3 && option
    previousMsje = {};
else
    previousMsje = get(handle,'String');
end
newMsje = {previousMsje{:} msje{:}}';
set(handle,'String',newMsje,'Value',numel(newMsje))
pause(0.1)
