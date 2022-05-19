k.mode <- function(data){
  k=density(data)
  mode=k[['x']][which.max(k[['y']])]
  return(mode)
}
find.spots <- function(data,box.size,low.lim,high.lim,fill.radius){
  peeks=splus2R::peaks(data,span = 2*box.size+1,strict = T)*t(splus2R::peaks(t(data),span = 2*box.size+1,strict = T))
  col.id=matrix(1:ncol(peeks),nrow = nrow(peeks),ncol = ncol(peeks),byrow = TRUE)
  row.id=matrix(1:nrow(peeks),nrow = nrow(peeks),ncol = ncol(peeks),byrow = FALSE)
  rows=c(row.id[peeks==1]); cols=c(col.id[peeks==1])
  vols=rep(0,times=sum(peeks)); for (i in 1:length(vols)){if(rows[i]>fill.radius & rows[i]<(nrow(data)-fill.radius) & cols[i]>fill.radius & cols[i]<(nrow(data)-fill.radius)){vols[i]=sum(data[(rows[i]-fill.radius):(rows[i]+fill.radius),(cols[i]-fill.radius):(cols[i]+fill.radius)])-median(data[(rows[i]-fill.radius):(rows[i]+fill.radius),(cols[i]-fill.radius):(cols[i]+fill.radius)])*(2*fill.radius+1)^2}}
  spots=list('x'=cols[(vols>=low.lim & vols<=high.lim)],'y'=((nrow(data)+1)-rows)[(vols>=low.lim & vols<=high.lim)],'row'=rows[(vols>=low.lim & vols<=high.lim)],'col'=cols[(vols>=low.lim & vols<=high.lim)],'vol'=vols[(vols>=low.lim & vols<=high.lim)])
  final=as.data.frame(cbind(spots[['x']],spots[['y']],spots[['row']],spots[['col']],spots[['vol']]));colnames(final)=col.names = c('x','y','row','col','vol')
  final=final[order(final$vol),]
  return(final)
}
classify.states <- function(data,signal.step=1000){
  fit=smooth.spline((1:length(data))[is.na(data)==F],na.omit(data),spar = 0.5)
  d1=predict(object = fit,x = 1:length(data),deriv = 1)
  swaps=c(1,d1[['x']][(abs(d1[['y']])>=1.0) & (splus2R::peaks(x = abs(d1[['y']]),span = 5))],max(d1[['x']]))
  pos.averages=rep(NA,times=length(data))
  states=rep(0,times=length(data))
  state.averages=rep(NA,times=length(data))
  for (i in 1:(length(swaps)-1)){
    pos.averages[swaps[i]:swaps[i+1]]=mean(na.omit(data[swaps[i]:swaps[i+1]]))
  }
  for (j in 0:50){
    states[(pos.averages>=(signal.step*j-signal.step/2)) & (pos.averages<=(signal.step*j+signal.step/2))]=j
  }
  for (k in 0:50){
    state.averages[states==k]=mean(na.omit(data[states==k]))
  }
  final=as.data.frame(cbind(states,state.averages))
  colnames(final)=c('id','avg')
  return(final)
}

id.spots <- function(path.to.file,file.name,time.step,spot.box=6,spot.radius=6,spot.min=300,spot.max=Inf){
  image.raw=suppressWarnings(tiff::readTIFF(paste0(path.to.file,file.name),all = TRUE,as.is=TRUE))
  pixel.size=dim(image.raw[[1]])
  frame.number=length(image.raw)
  image.data=array(NA,dim = c(pixel.size[1],pixel.size[2],frame.number)); for (i in 1:frame.number){image.data[,,i]=image.raw[[i]]}
  image.avg=apply(image.data,MARGIN = c(1,2),FUN = mean)
  #
  image(t(pracma::flipud(image.avg)),col=gray.colors(cumprod(pixel.size)),x = 1:pixel.size[2],y=1:pixel.size[1],axes=FALSE,xlab='',ylab='')
  check.1=readline(prompt = 'Image loading complete. Proceed with spot detection? (y/n):  '); if (check.1=='n'){stop('SCRIPT ABORTED BY USER')}
  #
  spots=find.spots(image.avg,spot.box,spot.min,spot.max,spot.radius)
  points(spots$x,spots$y,col='red',cex=1.2)
  show(paste0('Median particle intensity = ',median(spots$vol)))
  show(paste0('Number of particles identified = ',length(spots$x)))
  #
  check.2=readline(prompt = 'Is initial particle selection acceptable? (y/n):  '); while (check.2=='n'){
    spot.box=as.numeric(readline(prompt = 'New spot box size value (ENTER for default):  '))
    if (is.na(spot.box)==T){
      spot.box=6
    }
    spot.radius=as.numeric(readline(prompt = 'New spot radius value (ENTER for default):  '))
    if (is.na(spot.radius)==T){
      spot.radius=6
    }
    spot.min=as.numeric(readline(prompt = 'New spot minimum threshold (ENTER for default):  '))
    if (is.na(spot.min)==T){
      spot.min=500
    }
    spot.max=as.numeric(readline(prompt = 'New spot maximum threshold (ENTER for default):  '))
    if (is.na(spot.max)==T){
      spot.max=Inf
    }
    spots=find.spots(image.avg,spot.box,spot.min,spot.max,spot.radius)
    image(t(pracma::flipud(image.avg)),col=gray.colors(cumprod(pixel.size)),x = 1:pixel.size[2],y=1:pixel.size[1],axes=FALSE,xlab='',ylab='')
    points(spots$x,spots$y,col='red',cex=1.2)
    show(paste0('Median particle intensity = ',median(spots$vol)))
    show(paste0('Number of particles identified = ',length(spots$x)))
    check.2=readline(prompt = 'Is initial particle selection acceptable? (y/n/quit):  ')
    if (check.2=='quit'){
      stop()
    }
  }
  #
  particle.traces=matrix(0,nrow = frame.number,ncol = nrow(spots))
  particle.snaps=array(0,dim = c(2*spot.radius+1,2*spot.radius+1,20,nrow(spots)))
  snap.length=frame.number/20
  for (i in 1:nrow(spots)){
    particle.traces[,i]=apply(image.data[(spots$row[i]-spot.radius):(spots$row[i]+spot.radius),(spots$col[i]-spot.radius):(spots$col[i]+spot.radius),],MARGIN = c(3),FUN = sum)-apply(image.data[(spots$row[i]-spot.radius):(spots$row[i]+spot.radius),(spots$col[i]-spot.radius):(spots$col[i]+spot.radius),],MARGIN = c(3),FUN = median)*((2*spot.radius+1)^2)
    for (j in 1:20){
      particle.snaps[,,j,i]=apply(image.data[(spots$row[i]-spot.radius):(spots$row[i]+spot.radius),(spots$col[i]-spot.radius):(spots$col[i]+spot.radius),((j-1)*snap.length+1):(j*snap.length)],MARGIN = c(1,2),FUN = mean)
    }
  }
  particle.traces=particle.traces-k.mode(particle.traces)
  row.names(particle.traces)=time.step*(1:frame.number)
  #
  utils::write.table(spots,file = paste0(path.to.file,'initial_particle_summary.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
  utils::write.table(particle.traces,file = paste0(path.to.file,'initial_particle_traces.txt'),quote = FALSE,sep = '\t',row.names = TRUE,col.names = FALSE)
  save(list = c('image.avg','spots','particle.traces','particle.snaps','time.step','frame.number','pixel.size','spot.radius'),file = paste0(path.to.file,'Initial-Particle_Data.RData'))
  id.spots.output.files=list('composite image data'=image.avg,'initial particle summary'=spots,'particle traces'=particle.traces,'particle images data'=particle.snaps)
  return(id.spots.output.files)
}

refine.particles <- function(path.to.file,file.name='Initial-Particle_Data.RData',skip.manual='n',signal.step=1000,auto.filter='none'){
  load(file = paste0(path.to.file,file.name))
  #
  particle.trace.rolls=apply(particle.traces,MARGIN = c(2),FUN = data.table::frollmean,n=5,align='center')
  particles.to.keep=rep(FALSE,times=nrow(spots))
  residence.times=c('Particle','Start','Stop','Residence','State Signal')
  dwell.calls=c('Particle','State','Dwell')
  residence.calls=c('Particle','Bound','Residence')
  state.calls=matrix(NA,nrow = frame.number,ncol = nrow(spots))
  filtering=rep(F,times=nrow(spots))
  COUNTER=0
  for (i in 1:nrow(spots)){
    states=classify.states(particle.trace.rolls[,i],signal.step)
    state.calls[,i]=states$id
    dwell.calls=rbind(dwell.calls,cbind(rep(i,times=length(rle(state.calls[,i])[['values']])),rle(state.calls[,i])[['values']],rle(state.calls[,i])[['lengths']]*time.step))
    residence.calls=rbind(residence.calls,cbind(rep(i,times=length(rle(state.calls[,i]>=1)[['values']])),rle(state.calls[,i]>=1)[['values']],rle(state.calls[,i]>=1)[['lengths']]*time.step))
    if  (auto.filter=='unbound') {
      filtering[i]=(sum(state.calls[,i])==0)
    }
    if  (auto.filter=='all.stable') {
      filtering[i]=(length(table(state.calls[,i]))==1)
    }
  }
  #
  plot(frame.number*time.step,rowsum(as.matrix(state.calls>0))/nrow(spots),type='l',main='Photobleaching Check',xlab='Time (s)',ylab='Proportion of Bound Particles')
  temp.3=readline('Proceed with particle refinement? (y/n):  ')
  if (temp.3=='n'){
    stop()
  }
  #
  for (i in 1:nrow(spots)){
    if (skip.manual=='n' & filtering[i]==FALSE){
      par(fig=c(0,1,0.6,1))
      temp.1=matrix(particle.snaps[,,,i],nrow = (2*spot.radius+1),ncol = (2*spot.radius+1)*20)
      image(t(pracma::flipud(rbind(temp.1[,1:(ncol(temp.1)/2)],temp.1[,(ncol(temp.1)/2+1):ncol(temp.1)]))),col=gray.colors(length(temp.1)),axes=FALSE)
      par(fig=c(0,1,0,0.75),new = TRUE)
      plot((1:frame.number)*time.step,particle.trace.rolls[,i],type='l',main = 'Particle Trace',xlab='Time (s)',ylab = 'Signal')
      lines(((1:frame.number)*time.step)[is.na(states$avg)==F],states$avg[is.na(states$avg)==F],col='green',lwd=3)
      temp.2=readline('Should this particle be used for analysis? (y/n/quit):  ')
      if (temp.2=='y'){
        particles.to.keep[i]=TRUE
        COUNTER=COUNTER+1
        event.number=as.numeric(readline('How many binding events will you record for this particle?:  '))
        event.times=matrix(0,nrow=event.number,ncol = 5); event.times[,1]=rep(COUNTER,times=event.number)
        for (j in 1:event.number){
          show('Please click on the starting point then stopping point of a binding event.')
          event.times[j,2:3]=as.numeric(identify((1:frame.number)*time.step,particle.trace.rolls[,i],n=2)*time.step)
          event.times[j,4]=diff(event.times[j,2:3])
          event.times[j,5]=mean(na.omit(particle.trace.rolls[(event.times[j,2]/time.step):(event.times[j,3]/time.step),i]))
          arrows(x0 = event.times[j,2],x1 = event.times[j,3],y0=mean(na.omit(particle.trace.rolls[(event.times[j,2]/time.step):(event.times[j,3]/time.step),i])),y1=mean(na.omit(particle.trace.rolls[(event.times[j,2]/time.step):(event.times[j,3]/time.step),i])),angle = 90,code = 3,col = 'red')
        }
        residence.times=rbind(residence.times,event.times)
        dev.print(pdf,paste0(path.to.file,'Particle_Trace_',COUNTER,'.pdf'))
      }
      if (temp.2=='quit'){
        stop()
      }
    }
    if (skip.manual=='y' & i==1){
      temp.2='blank'
    }
  }
  residence.calls=as.data.frame(as.matrix(residence.calls[2:nrow(residence.calls),])); for (k in 1:3){residence.calls[,k]=as.numeric(residence.calls[,k])}; colnames(residence.calls)=c('Particle','Bound','Residence')
  dwell.calls=as.data.frame(as.matrix(dwell.calls[2:nrow(dwell.calls),])); for (k in 1:3){dwell.calls[,k]=as.numeric(dwell.calls[,k])}; colnames(dwell.calls)=c('Particle','State','Dwell')
  if (skip.manual=='n'){
    residence.data=as.data.frame(as.matrix(residence.times[2:nrow(residence.times),])); for (k in 1:5){residence.data[,k]=as.numeric(residence.data[,k])}; colnames(residence.data)=c('Particle','Start','Stop','Residence','State Signal')
    refined.particle.traces=particle.traces[,particles.to.keep]
    refined.state.calls=state.calls[,particles.to.keep]
    refined.particle.trace.rolls=particle.trace.rolls[,particles.to.keep]
    refined.spots=spots[particles.to.keep,]
    refined.particle.snaps=particle.snaps[,,,particles.to.keep]
  }
  #
  if (skip.manual=='n'){
    save(list = c('image.avg','residence.data','refined.particle.traces','refined.particle.trace.rolls','refined.spots','refined.particle.snaps','state.calls','residence.calls','dwell.calls'),file = paste0(path.to.file,'Refined-Particle_Data.RData'))
    utils::write.table(residence.data,file = paste0(path.to.file,'residence_data.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
    utils::write.table(refined.particle.traces,file = paste0(path.to.file,'selected_particle_traces.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
    utils::write.table(refined.particle.trace.rolls,file = paste0(path.to.file,'selected_particle_smoothed-traces.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
    utils::write.table(refined.spots,file = paste0(path.to.file,'selected_particle_summary.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
    return(list('image.avg'=image.avg,'residence.data'=residence.data,'refined.particle.traces'=refined.particle.traces,'refined.particle.trace.rolls'=refined.particle.trace.rolls,'refined.spots'=refined.spots,'refined.particle.snaps'=refined.particle.snaps,'state.calls'=state.calls,'residence.calls'=residence.calls,'dwell.calls'=dwell.calls))
  }
  if (skip.manual=='y'){
    save(list = c('image.avg','state.calls','residence.calls','dwell.calls'),file = paste0(path.to.file,'Refined-Particle_Data.RData'))
    return(list('image.avg'=image.avg,'state.calls'=state.calls,'residence.calls'=residence.calls,'dwell.calls'=dwell.calls))
  }
  utils::write.table(state.calls,file = paste0(path.to.file,'all-particle_state-calls.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
  utils::write.table(dwell.calls,file = paste0(path.to.file,'all-particle_dwell-calls.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
  utils::write.table(residence.calls,file = paste0(path.to.file,'all-particle_residence-calls.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
}

calc.kn1 <- function(path.to.file,file.name='Refined-Particle_Data.RData',use.auto.times='n',min.residence=0,max.residence=Inf){
  load(file = paste0(path.to.file,file.name))
  #
  if (use.auto.times=='n'){
    residence.times=residence.data$Residence
  }
  if (use.auto.times=='y'){
    residence.times=na.omit(residence.calls$Residence[residence.calls$Bound==T])
  }
  #
  residence.time.values=residence.times[((residence.times>=min.residence) & (residence.times<=max.residence))]
  residence.time.outliers=residence.times[((residence.times<min.residence) | (residence.times>max.residence))]
  #
  mod1=fitdistrplus::fitdist(data=as.numeric(residence.time.values),distr = 'gamma',method = c('mle'))
  mod1.shape=mod1$estimate[1]
  mod1.rate=mod1$estimate[2]
  #
  fit.x=seq(0,1.5*max(residence.time.values),0.01)
  mod1.fit=dgamma(fit.x,shape = mod1.shape,rate = mod1.rate)
  #
  par(fig =c(0,1,0.3,1))
  plot(fit.x,mod1.fit,type='l',xlab = 'Residence Time (s)',ylab = 'Probability Density',main = 'Distribution of Residence Times',col='blue')
  lines(density(residence.time.values),col='red')
  legend('topright',legend = c('Kernel Density','Gamma Fit','Data Points','Outliers'),fill = c('red','blue','black','grey'),col= c('red','blue','black','grey'))
  text(0.5*median(fit.x),max(mod1.fit),labels = paste0('Average Residence Time = ',signif(mean(residence.time.values),3),' (s)'),adj=c(0,1))
  text(0.5*median(fit.x),max(mod1.fit),labels = paste0('Dissociation Rate = ',signif(1/mean(residence.time.values),3),' (1/s)'),adj=c(0,3))
  #
  par(fig =c(0,1,0,0.4),new = TRUE)
  plot(residence.time.values,rep(0,times=length(residence.time.values)),type = 'p',col='black',main = NULL,xlim = range(fit.x),ann=FALSE,yaxt='n',xaxt='n')
  points(residence.time.outliers,rep(0,times=length(residence.time.outliers)),col='grey')
  abline(v=mean(residence.time.values),col='green',lwd=3)
  #
  show(paste0('Average Residence Time = ',signif(mean(residence.time.values),3),' s'))
  show(paste0('Dissociation Rate = ',signif(1/mean(residence.time.values),3),' 1/s'))
  #
  check.1=readline('Do you wish to save this analysis? (y/n):   ')
  if (check.1=='y'){
    dev.print(pdf,paste0(path.to.file,'Dissociation-Rate_Analysis.pdf'))
  }
}






