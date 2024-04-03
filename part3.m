clc, clear, close all

load('zpw.mat')
load('OFDM_PILOT.mat')
load('ofdm_map.mat')

k = 2048;
L = 199;
kp = 512;
null_sub = 112;

idx = find(ofdm_map ==1);


for nn = 1:kp
    for n = 1:L+1
        V(nn,n) = exp((-1i*2*pi*(n-1)*(idx(nn)-1))/k);
    end
end

D = OFDM_PILOT(ofdm_map==1);

D = diag(D);
hls = (1/kp)*V'*D'*Z_pw;

for nn = 1:kp
    for n = 1:L+1
        V2(nn,n) = exp((-1i*2*pi*(n-1)*(nn-1))/k);
    end
end

Hw = V2*hls;

Z_pw_null = Z_pw(ofdm_map ==0);

noiseVarience = (1/null_sub);