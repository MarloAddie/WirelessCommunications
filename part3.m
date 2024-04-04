clc, clear, close all

load('zpw.mat')
load('OFDM_PILOT.mat')
load('ofdm_map.mat')
load('z_w.mat')

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

for nn = 1:k
    for n = 1:L+1
        V2(nn,n) = exp((-1i*2*pi*(n-1)*(nn-1))/k);
    end
end

Hw = V2*hls;

Z_pw_null = z_w((ofdm_map == 0));

noiseVar = (1/null_sub)*sum(abs(Z_pw_null).^2);

sd = find(ofdm_map == 2);

zd = z_w(sd-1);

HD = Hw(sd-1);

X1 = 1/sqrt(2) + (1/sqrt(2))*1i;
X2 = -1/sqrt(2) + (1/sqrt(2))*1i;
X3 = 1/sqrt(2) - (1/sqrt(2))*1i;
X4 = -1/sqrt(2) - (1/sqrt(2))*1i;

LR = [];
for pl = 1:length(HD)
    
    Lb1 = log10((exp(-norm(zd(pl)-HD(pl)*X1)^2/noiseVar)+norm(-abs(zd(pl)-HD(pl)*X3)^2/noiseVar))...
        /(exp(-norm(zd(pl)-HD(pl)*X2)^2/noiseVar)+exp(-norm(zd(pl)-HD(pl)*X4)^2/noiseVar)));
    
    Lb2 = log10((exp(-norm(zd(pl)-HD(pl)*X1)^2/noiseVar)+exp(-norm(zd(pl)-HD(pl)*X2)^2/noiseVar))...
        /(exp(-norm(zd(pl)-HD(pl)*X3)^2/noiseVar)+exp(-norm(zd(pl)-HD(pl)*X4)^2/noiseVar)));
    
    LR1 = cat(1, Lb1, Lb2);
    
    LR = cat(1, LR, LR1);
    
end





