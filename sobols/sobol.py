from networkx import connected_watts_strogatz_graph as nxs
from networkx import powerlaw_cluster_graph as nxl
from itertools import starmap, product, combinations
from time import time
from matplotlib import pyplot
from multiprocessing import Pool
import os
from numpy import sqrt, mean, array, matrix, ndarray, argsort, tensordot, sort, ones, zeros
from numpy.linalg import norm
from numpy.random import choice, rand

states = [ 'S' , 'A' , 'I' , 'R' ]
colors = { 'S' : 'tab:blue' , 'A' : 'tab:orange' , 'I' : 'tab:green' , 'R' : 'tab:red' }
back_reactions = { 'S' : [ 'A' , 'I' ] , 'A' : [] , 'I' : [] , 'R' : [] }
forw_reactions = { 'S' : [] ,'A' : [ 'S' ] , 'I' : [ 'S' ] , 'R' : [] }
reaction = [ 'alpha' , 'betaA' , 'betaI' , 'gamma' , 'detaA' , 'detaI' , 'varnu' ]
interval = { 'alpha' : ( 0.001 , 0.9 ) , 'betaA' : ( 0.5 , 0.9 ) , 'betaI' : ( 0.5 , 0.9 ) , 'gamma' : ( 0.001 , 0.05 ) , 'detaA' : ( 0.1 , 0.51 ) , 'detaI' : ( 0.1 , 0.51 ) , 'varnu' : ( 0 , 0.02 ) }
befor = { 'alpha' : 'A' , 'varnu' : 'S' , 'gamma' : 'R' , 'detaA' : 'A' , 'detaI' : 'I' }
befor |= { 'betaA' : ( 'A' , 'S' ) , 'betaI' : ( 'I' , 'S' ) }
after = { 'alpha' : 'I' , 'varnu' : 'R' , 'gamma' : 'S' , 'detaA' : 'R' , 'detaI' : 'R' }
after |= { 'betaA' : 'A' , 'betaI' : 'A' }

def update( k , state , elem , x , contact , num ):
	for h in contact[ k ]:
		if ( x[ h ] in back_reactions[ x[ k ] ] ):
			num[ ( x[ h ] , x[ k ] ) ] -= 1
			elem[ ( x[ h ] , x[ k ] ) ].remove( ( h , k ) )
		if ( x[ h ] in forw_reactions[ x[ k ] ] ):
			num[ ( x[ k ] , x[ h ] ) ] -= 1
			elem[ ( x[ k ] , x[ h ] ) ].remove( ( k , h ) )
		if ( x[ h ] in back_reactions[ state ] ):
			num[ ( x[ h ] , state ) ] += 1
			elem[ ( x[ h ] , state ) ] += [ ( h , k ) ]
		if ( x[ h ] in forw_reactions[ state ] ):
			num[ ( state , x[ h ] ) ] += 1
			elem[ ( state , x[ h ] ) ] += [ ( k , h ) ]
	num[ x[ k ] ] -= 1
	num[ state ] += 1
	elem[ x[ k ] ].remove( k )
	elem[ state ] += [ k ]
	x[ k ] = state

def init( citizens , net ):
	dgs = zeros( citizens )
	max = 1
	arg = []
	for ( k , h ) in net:
		dgs[ k ] += 1
		dgs[ h ] += 1
		if ( dgs[ k ] == max ):
			arg += [ k ]
			max = dgs[ k ]
		if ( dgs[ k ] > max ):
			arg = [ k ]
			max = dgs[ k ]
		if ( dgs[ h ] == max ):
			arg += [ h ]
			max = dgs[ h ]
		if ( dgs[ h ] > max ):
			arg = [ h ]
			max = dgs[ h ]
	temp = array( citizens * [ 'S' ] )
	for k in arg: temp[ k ] = 'A'
	return temp

def stochastic( net , Xparam , days , citizens ):
	x = init( citizens , net )
	num = {}
	elem = {}
	for rea in reaction:
		num[ befor[ rea ] ] = 0
		elem[ befor[ rea ] ] = []
	contact = []
	for k in range( citizens ):
		contact += [ [] ]
		num[ x[ k ] ] += 1
		elem[ x[ k ] ] += [ k ]
	for ( k , h ) in net:
		contact[ k ] += [ h ]
		contact[ h ] += [ k ]
		if x[ h ] in back_reactions[ x[ k ] ]:
			elem[ ( x[ h ] , x[ k ] ) ] += [ ( h , k ) ]
			num[ ( x[ h ] , x[ k ] ) ] += 1
		if x[ k ] in back_reactions[ x[ h ] ]:
			elem[ ( x[ k ] , x[ h ] ) ] += [ ( k , h ) ]
			num[ ( x[ k ] , x[ h ] ) ] += 1
	y = { 'time' : [ 0 ] }
	for state in states: y[ state ] = [ num[ state ] ]
	t = 0
	while ( t < days ):
		randomic = rand()
		prob = {}
		rate = 0
		for rea in reaction:
			prob[ rea ] = num[ befor[ rea ] ] * Xparam[ rea ]
			rate += prob[ rea ]
		dt = 1 / rate
		t += dt
		for rea in reaction:
			prob[ rea ] *= dt
		pro = 0
		count = 0
		while ( randomic > pro + prob[ reaction[ count ] ] ):
			pro += prob[ reaction[ count ] ]
			count += 1
		rea = reaction[ count ]
		k = int( float( num[ befor[ rea ] ] ) * ( randomic - pro ) / prob[ rea ] )
		if ( befor[ rea ] in states ):
			k = elem[ befor[ rea ] ][ k ]
		else:
			h , k = elem[ befor[ rea ] ][ k ]
		update( k , after[ rea ] , elem , x , contact , num )
		for state in states:
			y[ state ] += [ num[ state ] ]
		y[ 'time' ] += [ t ]
		if ( num[ 'A' ] == num[ 'I' ] == 0 ):
			t = days
			y[ 'time' ][ - 1 ] = days
	return y

def fold( path ):
	if not os.path.exists( path ):
		try: os.mkdir( path )
		except: return

def form( some ):
	if ( type( some ) == list ):
		TEMP = '['
		for temp in some:
			TEMP += form( temp )
		return TEMP + '],'
	if ( type( some ) == matrix ) or ( type( some ) == ndarray ):
		return form( list( array( some ) ) )
	if type( some ) == tuple:
		TEMP = '('
		for temp in some:
			TEMP += form( temp )
		return TEMP + '),'
	return str( some ) + ','

def save( some , path ):
	f = open( path , 'w' )
	f.write( form( some ) )
	f.close()

def load( path ):
	try: f = open( path , 'r' )
	except: return []
	TEMP = [[]]
	temp = ''
	while( True ):
		a = f.read( 1 )
		if ( a == '' ):
			f.close()
			if( len( TEMP ) == len( TEMP[ 0 ] ) == 1 ): return TEMP[0][0]
			else: raise Exception()
		elif ( a == '[' ): TEMP += [[]]
		elif ( a == '(' ): TEMP += [[]]
		elif ( a == ']' ): temp = [] + TEMP.pop( -1 )
		elif ( a == ')' ):
			temp = TEMP.pop( -1 )
			if len( temp ) == 2:
				temp = ( temp[ 0 ] , temp[ 1 ] )
			else:
				f.close()
				f = open( path , 'r' )
				print( f.read() )
				print( TEMP )
				f.close()
				raise Exception()
		elif ( a == ',' ):
			try:
				temp = float( temp )
				if ( temp == int( temp ) ): TEMP[ -1 ] += [ int( temp ) ]
				else: TEMP[ -1 ] += [ temp ]
			except:
				TEMP[ -1 ] += [ temp ]
			temp = ''
		else: temp += a

def draw( y , ax , citizens , lw ):
	ax[ 1 , 0 ].set_xlabel( 'days' )
	ax[ 1 , 1 ].set_xlabel( 'days' )
	ax[ 0 , 0 ].set_ylabel( 'individuals' )
	ax[ 1 , 0 ].set_ylabel( 'individuals' )
	for sta in [ 0 , 1 , 2 , 3 ]:
		ax[ sta % 2 , int( sta / 2 ) ].plot( y[ 'time' ] , array( y[ states[ sta ] ] ) / citizens , lw = lw , c = colors[ states[ sta ] ] )
		ax[ sta % 2 , int( sta / 2 ) ].legend( [ states[ sta ] ] , loc = 'upper right' )

class mydict( dict ):
	def __mul__( self , x ): return mydict( self | x )

class Admin:
	def __init__( self , net , citizens , days ):
		fold( '../Data' )
		self.folder = '../Data/' + net + '_' + str( citizens )
		fold( self.folder )
		self.topology = net[ 0 : 3 ]
		self.deg = int( net[ 3 ] )
		self.citizens = citizens
		self.days = days
		self.param = {}
		for par in reaction: self.param[ par ] = load( self.folder + '/' + par )
		self.theta = load( self.folder + '/theta' )

	def newpar( self , par ):
		if ( par == 'theta' ):
			self.theta += [ rand() ]
			save( self.theta , self.folder + '/theta' )
		else:
			self.param[ par ] += [ interval[ par ][ 0 ] + rand() * ( interval[ par ][ 1 ] - interval[ par ][ 0 ] ) ]
			save( self.param[ par ] , self.folder + '/' + par )

	def newnet( self , k ):
		while ( len( self.theta ) <= k ): self.newpar( 'theta' )
		fold( self.folder + '/' + str( k ) )
		see = load( self.folder + '/' + str( k ) + '_see' )
		seed = len( see )
		save( see + [[]] , self.folder + '/' + str( k ) + '_see' )
		if self.topology == 'nxs': net = list( nxs( self.citizens , self.deg , self.theta[ k ] ).edges() )
		if self.topology == 'nxl': net = list( nxl( self.citizens , self.deg , self.theta[ k ] ).edges() )
		save( net , self.folder + '/' + str( k ) + '/' + str( seed ) + '_net' )
		fold( self.folder + '/' + str( k ) + '/' + str( seed ) )

	def newsim( self , param , net ):
		Xparam = {} | param
		par_id = ''
		for par in reaction:
			while ( len( self.param[ par ] ) <= Xparam[ par ] ): self.newpar( par )
			par_id += str( Xparam[ par ] ) + '_'
			Xparam[ par ] = self.param[ par ][ Xparam[ par ] ]
		temp = stochastic( net , Xparam , self.days , self.citizens )
		infe = array( temp[ 'A' ] ) + array( temp[ 'I' ] )
		if ( 2 * temp[ 'time' ][ infe.argmax() ] > self.days ): print( 'Warning: peak not found: days =' , self.days )
		if ( 4 * temp[ 'time' ][ infe.argmax() ] > self.days ): self.days *= 2
		self.days += temp[ 'time' ][ infe.argmax() ] - self.days / 8
		return infe.max()
		
	def mytest( self ):
		Xparam = { 'theta' : choice( range( len( self.theta ) ) ) }
		for par in reaction: Xparam[ par ] = choice( range( len( self.param[ par ] ) ) )
		par_id = ''
		for par in reaction:
			par_id += str( Xparam[ par ] ) + '_'
			Xparam[ par ] = self.param[ par ][ Xparam[ par ] ]
		see = load( self.folder + '/' + str( Xparam[ 'theta' ] ) + '_see' )
		seed = choice( range( len( see ) ) )
		net = load( self.folder + '/' + str( Xparam[ 'theta' ] ) + '/' + str( seed ) + '_net' )
		fig , ax = pyplot.subplots( 2 , 2 , squeeze = False )
		draw( stochastic( net , Xparam , self.days , self.citizens ) , ax , self.citizens , 1 )
		fig.tight_layout()
		file = self.folder + '/' + str( Xparam[ 'theta' ] ) + '/' + str( seed )
		fold( file )
		file += '/' + par_id[ 0 : -1 ]
		fold( file )
		file += '/' + str( self.days )
		fold( file )
		pyplot.savefig( file + '/' + str( len( os.listdir( file ) ) ) + '_plot.png' )
		pyplot.close( fig )

	def newmat( self , param , fill ):
		Xparam = {} | param
		par_id = ''
		for par in reaction: par_id += str( Xparam[ par ] ) + '_'
		if ( len( fill ) == 3 ):
			try:
				mat = array( load( self.folder + '/' + par_id + 'mat' ) )
				fill0 = sort( choice( mat.shape[ 0 ] , fill[ 0 ] , replace = False ) )
				fill1 = sort( choice( mat.shape[ 1 ] , fill[ 1 ] , replace = False ) )
				fill2 = sort( choice( mat.shape[ 2 ] , fill[ 2 ] , replace = False ) )
				mat = mat[ fill0 ]
				mat = array( list( map( lambda x : matrix( x )[ fill1 , : ][ : , fill2 ] , list( mat ) ) ) )
			except:
				mat = list( starmap( self.newmat , map( lambda x : ( Xparam | { 'theta' : x } , [ fill[ 1 ] , fill[ 2 ] ] ) , range( fill[ 0 ] ) ) ) )
				mat = array( mat )
				save( mat , self.folder + '/' + par_id + 'mat' )
			return mat[ Xparam[ 'theta' ] ]
		if ( len( fill ) == 2 ):
			try:
				mat = matrix( load( self.folder + '/' + str( Xparam[ 'theta' ] ) + '/' + par_id + 'mat' ) )
				fill0 = sort( choice( mat.shape[ 0 ] , fill[ 0 ] , replace = False ) )
				fill1 = sort( choice( mat.shape[ 1 ] , fill[ 1 ] , replace = False ) )
				mat = mat[ fill0 , : ][ : , fill1 ]
			except:
				mat = list( starmap( self.newmat , map( lambda x : ( Xparam | { 'lamb0' : x } , [ fill[ 1 ] ] ) , range( fill[ 0 ] ) ) ) )
				mat = sort( array( mat ) , axis = 0 )
				save( mat , self.folder + '/' + str( Xparam[ 'theta' ] ) + '/' + par_id + 'mat' )
			return mat[ Xparam[ 'lamb0' ] ]
		if ( len( fill ) == 1 ):
			try:
				mat = array( load( self.folder + '/' + str( Xparam[ 'theta' ] ) + '/' + str( Xparam[ 'lamb0' ] ) + '/' + par_id + 'mat' ) )
				fill0 = sort( choice( mat.shape[ 0 ] , fill[ 0 ] , replace = False ) )
				mat = mat[ fill0 ]
			except:
				net = load( self.folder + '/' + str( Xparam[ 'theta' ] ) + '/' + str( Xparam[ 'lamb0' ] ) + '_net' )
				mat = list( starmap( self.newsim , map( lambda x : ( Xparam | { 'lamb1' : x } , net ) , range( fill[ 0 ] ) ) ) )
				mat = sort( array( mat ) )
				save( mat , self.folder + '/' + str( Xparam[ 'theta' ] ) + '/' + str( Xparam[ 'lamb0' ] ) + '/' + par_id + 'mat' )
			return mat[ Xparam[ 'lamb1' ] ]

	def mysobo( self , param , fill ):
		for k in range( fill[ 0 ] ):
			while ( len( load( self.folder + '/' + str( k ) + '_see' ) ) < fill[ 1 ] ): self.newnet( k )
		Xparam = mydict( param )
		for par in reaction: Xparam = tensordot( Xparam , array( list( map( lambda x : { par : x } , param[ par ] ) ) ) , axes = 0 )
		shape = Xparam.shape
		with Pool() as pool: Xparam = list( pool.starmap( self.newmat , map( lambda x : ( x , fill ) , Xparam.reshape( Xparam.size ) ) , chunksize = 1 ) )
		Xparam = array( list( Xparam ) ).reshape( shape + ( len( param[ 'theta' ] ) , len( param[ 'lamb0' ] ) , len( param[ 'lamb1' ] ) ) )
		pick = []
		for r in range( 11 ): pick += list( combinations( range( 10 ) , r ) )
		freeze = list( map( lambda x : mean( mean( Xparam , axis = x )**2 ) , pick ) )
		freeze.reverse()
		X = {}
		Y = zeros( 10 )
		for p , v in zip( pick , freeze ):
			temp = v
			for k in range( len( p ) ):
				for q in combinations( p , k ): temp -= X[ q ]
			X[ p ] = temp
			for a in p: Y[ a ] += temp
		return Y

if ( __name__ == '__main__' ):
	M = 3
	param = {}
	for par in reaction + [ 'theta' , 'lamb0' , 'lamb1' ]: param |= { par : list( range( M ) ) }
	ad = Admin( 'nxl3' , 8192 , 16 )
	print( ad.mysobo( param , [ M , M , M ] ) )
