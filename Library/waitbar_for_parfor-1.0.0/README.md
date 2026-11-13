This script generates a waitbar or a wait message including (1) the percentage finished; (2) the time used; and (3) the estimated remaining time for parfor or other types of loop. The script uses "parallel.pool.DataQueue" to collect information from workers.

Example.m provides an example. Basically, we need to construct the object right before the loop, send a message to the object within the loop, and destroy the object after the loop.

%Example:  
N = 100000;  
WaitMessage = parfor\_wait(N, 'Waitbar', true);  
parfor i = 1: N  
 WaitMessage.Send;  
 pause(0.002);  
end  
WaitMessage.Destroy