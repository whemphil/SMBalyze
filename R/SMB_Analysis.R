id.spots <- function(path.to.file,file.name,time.step,spot.box=6,spot.radius=6,spot.min=300,spot.max=Inf,spot.picking='composite',r.pick=2.9,min.pick=5){
  find.spots <- function(data,box.size,low.lim,high.lim,fill.radius,spot.picking,r.pick,min.pick){
    if (spot.picking=='composite'){
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
    if (spot.picking=='all.frames'){
      spots=list(NULL)
      peeks=splus2R::peaks(data,span = 2*box.size+1,strict = T)*aperm(splus2R::peaks(aperm(data,c(2,1,3)),span = 2*box.size+1,strict = T),c(2,1,3))
      col.id=aperm(array(1:dim(peeks)[2],dim = dim(peeks)[c(2,1,3)]),c(2,1,3))
      row.id=array(1:dim(peeks)[1],dim = dim(peeks))
      layer.id=array(rep(1:dim(peeks)[3],each=dim(peeks)[1]*dim(peeks)[2]),dim = dim(peeks))
      rows=c(row.id[peeks==1]); cols=c(col.id[peeks==1]); layers=c(layer.id[peeks==1])
      vols=rep(0,times=sum(peeks)); for (j in 1:length(vols)){if(rows[j]>fill.radius & rows[j]<(dim(data)[1]-fill.radius) & cols[j]>fill.radius & cols[j]<(dim(data)[1]-fill.radius)){vols[j]=sum(data[(rows[j]-fill.radius):(rows[j]+fill.radius),(cols[j]-fill.radius):(cols[j]+fill.radius),layers[j]])-median(data[(rows[j]-fill.radius):(rows[j]+fill.radius),(cols[j]-fill.radius):(cols[j]+fill.radius),layers[j]])*(2*fill.radius+1)^2}}
      spots=data.frame('frame'=layers[(vols>=low.lim & vols<=high.lim)],'x'=cols[(vols>=low.lim & vols<=high.lim)],'y'=((dim(data)[1]+1)-rows)[(vols>=low.lim & vols<=high.lim)],'row'=rows[(vols>=low.lim & vols<=high.lim)],'col'=cols[(vols>=low.lim & vols<=high.lim)],'vol'=vols[(vols>=low.lim & vols<=high.lim)])
      clusters=dbscan::dbscan(spots[,2:3],eps = r.pick,minPts = min.pick)[['cluster']]
      particles=aggregate(x=as.data.frame(cbind(spots[,2:6],'cluster'=clusters)),by = cluster,FUN = mean)
      particles=particles[particles$cluster!=0,]; particles=particles[order(particles$vol),]; particles$col=round(particles$col); particles$row=round(particles$row)
      return(particles)
    }
  }
  image.raw=suppressWarnings(tiff::readTIFF(paste0(path.to.file,file.name),all = TRUE,as.is=TRUE))
  pixel.size=dim(image.raw[[1]])
  frame.number=length(image.raw)
  image.data=array(NA,dim = c(pixel.size[1],pixel.size[2],frame.number)); for (i in 1:frame.number){image.data[,,i]=image.raw[[i]]}
  image.avg=apply(image.data,MARGIN = c(1,2),FUN = mean)
  #
  suppressWarnings(image(t(pracma::flipud(image.avg)),col=gray.colors(cumprod(pixel.size)),x = 1:pixel.size[2],y=1:pixel.size[1],axes=FALSE,xlab='',ylab=''))
  check.1=readline(prompt = 'Image loading complete. Proceed with spot detection? (y/n):  '); if (check.1=='n'){stop('SCRIPT ABORTED BY USER')}
  #
  if (spot.picking=='composite'){
    spots=find.spots(image.avg,spot.box,spot.min,spot.max,spot.radius,spot.picking,r.pick,min.pick)
  }
  if (spot.picking=='all.frames'){
    spots=find.spots(image.data,spot.box,spot.min,spot.max,spot.radius,spot.picking,r.pick,min.pick)
  }
  points(spots$x,spots$y,col='red',cex=1.2)
  show(paste0('Median particle intensity = ',median(spots$vol)))
  show(paste0('Number of particles identified = ',length(spots$x)))
  #
  check.2=readline(prompt = 'Is initial particle selection acceptable? (y/n):  '); while (check.2=='n'){
    spot.box=as.numeric(readline(prompt = 'New spot selection box size (ENTER for default):  '))
    if (is.na(spot.box)==T){
      spot.box=6
    }
    spot.radius=as.numeric(readline(prompt = 'New spot integration radius (ENTER for default):  '))
    if (is.na(spot.radius)==T){
      spot.radius=6
    }
    spot.min=as.numeric(readline(prompt = 'New spot minimum signal (ENTER for default):  '))
    if (is.na(spot.min)==T){
      spot.min=300
    }
    spot.max=as.numeric(readline(prompt = 'New spot maximum signal (ENTER for default):  '))
    if (is.na(spot.max)==T){
      spot.max=Inf
    }
    if (spot.picking=='all.frames'){
      r.pick=as.numeric(readline(prompt = 'New spot merging radius (ENTER for default):  '))
      if (is.na(r.pick)==T){
        r.pick=2.9
      }
      min.pick=as.numeric(readline(prompt = 'New spot number minimum (ENTER for default):  '))
      if (is.na(min.pick)==T){
        min.pick=5
      }
    }
    if (spot.picking=='composite'){
      spots=find.spots(image.avg,spot.box,spot.min,spot.max,spot.radius,spot.picking,r.pick,min.pick)
    }
    if (spot.picking=='all.frames'){
      spots=find.spots(image.data,spot.box,spot.min,spot.max,spot.radius,spot.picking,r.pick,min.pick)
    }
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
  row.names(particle.traces)=time.step*(1:frame.number)
  #
  utils::write.table(spots,file = paste0(path.to.file,'initial_particle_summary.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
  utils::write.table(particle.traces,file = paste0(path.to.file,'initial_particle_traces.txt'),quote = FALSE,sep = '\t',row.names = TRUE,col.names = FALSE)
  save(list = c('image.avg','spots','particle.traces','particle.snaps','time.step','frame.number','pixel.size','spot.radius'),file = paste0(path.to.file,'Initial-Particle_Data.RData'))
  id.spots.output.files=list('composite image data'=image.avg,'initial particle summary'=spots,'particle traces'=particle.traces,'particle images data'=particle.snaps)
  return(id.spots.output.files)
}

refine.particles <- function(path.to.file,file.name='Initial-Particle_Data.RData',skip.manual='n',signal.step=NULL,auto.filter='none',classification.strategy='classic',background.subtraction='lower.quartile'){
  k.mode <- function(data){
    k=density(data)
    mode=k[['x']][which.max(k[['y']])]
    return(mode)
  }
  classify.states <- function(data,signal.step=1000,fun='classic'){
    if (fun=='classic'){
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
    if (fun=='basic'){
      pos.averages=matrix(NA,ncol = ncol(data),nrow = nrow(data))
      for (n in 1:ncol(data)){
        fit=smooth.spline((1:length(data[,n]))[is.na(data[,n])==F],na.omit(data[,n]),spar = 0.5)
        d1=predict(object = fit,x = 1:length(data[,n]),deriv = 1)
        swaps=c(1,d1[['x']][(abs(d1[['y']])>=1.0) & (splus2R::peaks(x = abs(d1[['y']]),span = 5))],max(d1[['x']]))
        for (i in 1:(length(swaps)-1)){
          pos.averages[swaps[i]:swaps[i+1],n]=mean(na.omit(data[swaps[i]:swaps[i+1],n]))
        }
      }
      final=pos.averages
      return(final)
    }
    if (fun=='uni.dbscan'){
      pos.averages=matrix(NA,ncol = ncol(data),nrow = nrow(data))
      for (n in 1:ncol(data)){
        fit=smooth.spline((1:length(data[,n]))[is.na(data[,n])==F],na.omit(data[,n]),spar = 0.5)
        d1=predict(object = fit,x = 1:length(data[,n]),deriv = 1)
        swaps=c(1,d1[['x']][(abs(d1[['y']])>=1.0) & (splus2R::peaks(x = abs(d1[['y']]),span = 5))],max(d1[['x']]))
        for (i in 1:(length(swaps)-1)){
          pos.averages[swaps[i]:swaps[i+1],n]=mean(na.omit(data[swaps[i]:swaps[i+1],n]))
        }
      }
      uni.dbscan <- function(datter,eps){
        temp.0=order(c(datter))
        temp.1=datter[temp.0]
        temp.2=c(abs(temp.1[1:(length(temp.1)-1)]-temp.1[2:length(temp.1)]),0)>eps
        temp.3=c(1,(1:length(temp.2))[temp.2],length(temp.2))
        temp.4=rep(NA,times=length(temp.2))
        for (i in 0:(length(temp.3)-2)){
          temp.4[(temp.3[i+1]):(temp.3[i+2])]=i
        }
        temp.5=matrix(temp.4[order(temp.0)],nrow = nrow(datter),ncol = ncol(datter))
        return(temp.5)
      }
      classes=uni.dbscan(pos.averages,eps = 1*sd(na.omit(data-pos.averages)))
      values=matrix(NA,nrow = nrow(classes),ncol = ncol(classes))
      for (k in 0:max(classes)){
        values[classes==k]=mean(na.omit(data[classes==k]))
      }
      final=list('id'=classes,'avg'=values)
      return(final)
    }
    if (fun=='var.shift'){
      pos.averages=rep(NA,times=length(data))
      fit=smooth.spline((1:length(data))[is.na(data)==F],na.omit(data),spar = 0.5)
      d1=predict(object = fit,x = 1:length(data),deriv = 1)
      swaps=c(1,d1[['x']][(abs(d1[['y']])>=1.0) & (splus2R::peaks(x = abs(d1[['y']]),span = 5))],max(d1[['x']]))
      for (i in 1:(length(swaps)-1)){
        pos.averages[swaps[i]:swaps[i+1]]=mean(na.omit(data[swaps[i]:swaps[i+1]]))
      }
      var.shift <- function(datter,eps){
        temp.0=order(datter)
        temp.1=datter[temp.0]
        temp.2=c(abs(temp.1[1:(length(temp.1)-1)]-temp.1[2:length(temp.1)]),0)>eps
        temp.3=c(1,(1:length(temp.2))[temp.2],length(temp.2))
        temp.4=rep(NA,times=length(temp.2))
        for (i in 1:(length(temp.3)-1)){
          temp.4[(temp.3[i]):(temp.3[i+1])]=i
        }
        temp.5=temp.4[order(temp.0)]
        return(temp.5)
      }
      classes=var.shift(pos.averages,eps = 1.7*sd(na.omit(data-pos.averages)))
      values=rep(NA,times=length(classes))
      for (k in 0:max(classes)){
        values[classes==k]=mean(na.omit(data[classes==k]))
      }
      classes[abs(values)<=2*sd(na.omit(data-pos.averages))]=0
      final=data.frame('id'=classes,'avg'=values)
      return(final)
    }
  }
  load(file = paste0(path.to.file,file.name))
  #
  particle.trace.rolls=apply(particle.traces,MARGIN = c(2),FUN = data.table::frollmean,n=5,align='center')
  if (background.subtraction=='k.mode'){
    particle.trace.rolls=particle.trace.rolls-k.mode(na.omit(particle.trace.rolls))
  }
  if (background.subtraction=='basal.states'){
    particle.trace.rolls=particle.trace.rolls-median(na.omit(apply(classify.states(particle.trace.rolls,fun = 'basic'),MARGIN = c(2),FUN = min)))
  }
  if (background.subtraction=='lower.quartile'){
    particle.trace.rolls=particle.trace.rolls-(c(na.omit(classify.states(particle.trace.rolls,fun = 'basic')))[order(c(na.omit(classify.states(particle.trace.rolls,fun = 'basic'))))])[round(0.25*length(c(na.omit(classify.states(particle.trace.rolls,fun = 'basic')))))]
  }
  particles.to.keep=rep(FALSE,times=nrow(spots))
  residence.times=c('Particle','Start','Stop','Residence','State Signal')
  dwell.calls=c('Particle','State','Dwell')
  residence.calls=c('Particle','Bound','Residence')
  if (classification.strategy=='classic' | classification.strategy=='var.shift'){
    state.calls=matrix(NA,nrow = frame.number,ncol = nrow(spots))
    if (is.null(signal.step)==TRUE){
      signal.step=2.5*sd(na.omit(particle.trace.rolls-classify.states(particle.trace.rolls,fun = 'basic')))
    }
  }
  if (classification.strategy=='uni.dbscan'){
    classifications=classify.states(particle.trace.rolls,fun = 'uni.dbscan')
    state.calls=classifications[['id']]
  }
  filtering=rep(F,times=nrow(spots))
  COUNTER=0
  for (i in 1:nrow(spots)){
    if (classification.strategy=='classic' | classification.strategy=='var.shift'){
      states=classify.states(particle.trace.rolls[,i],signal.step,fun = classification.strategy)
      state.calls[,i]=states$id
    }
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
  plot((1:frame.number)*time.step,apply(state.calls>0,MARGIN = c(1),FUN = sum)/nrow(spots),type='l',main='Photobleaching Check',xlab='Time (s)',ylab='Proportion of Bound Particles')
  temp.3=readline('Proceed with particle refinement? (y/n):  ')
  if (temp.3=='n'){
    stop()
  }
  #
  for (i in 1:nrow(spots)){
    if (skip.manual=='n' & filtering[i]==FALSE){
      if (classification.strategy=='classic' | classification.strategy=='var.shift'){
        states=classify.states(particle.trace.rolls[,i],signal.step,fun = classification.strategy)
      }
      if (classification.strategy=='uni.dbscan'){
        states=data.frame('id'=classifications[['id']][,i],'avg'=classifications[['avg']][,i])
      }
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






