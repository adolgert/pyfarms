

def make_map():
    try:
        import mapnik
    except:
        import mapnik
    m=mapnik.Map(600, 300)
    m.background=mapnik.Color("steelblue")

    s = mapnik.Style() # style object to hold rules
    r = mapnik.Rule() # rule object to hold symbolizers
    # to fill a polygon we create a PolygonSymbolizer
    polygon_symbolizer = mapnik.PolygonSymbolizer(mapnik.Color('#f2eff9'))
    r.symbols.append(polygon_symbolizer) # add the symbolizer to the rule object
    # to add outlines to a polygon we create a LineSymbolizer
    line_symbolizer = mapnik.LineSymbolizer(mapnik.Color('rgb(50%,50%,50%)'),0.1)
    r.symbols.append(line_symbolizer) # add the symbolizer to the rule object
    s.rules.append(r) # now add the rule to the style and we're done

    m.append_style('My Style',s) # Styles are given names only as they are applied to the map

    # South Carolina 400 Scale Grid. http://gis.sc.gov/data.html
    ds = mapnik.Shapefile(file='400/400.shp')
    print(ds.envelope())
    # Box2d(1290000.0,70000.0,2750000.0,1240000.0)

    layer = mapnik.Layer('SouthCarolina') # new layer called 'world' (we could name it anything)
    layer.srs='+proj=lcc +lat_1=32.5 +lat_2=34.83333333333334 +lat_0=31.83333333333333 +lon_0=-81 +x_0=609600 +y_0=0 +datum=NAD83 +units=ft +no_defs '
    # note: layer.srs will default to '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    layer.datasource = ds

    layer.styles.append('My Style')
    m.layers.append(layer)
    m.zoom_all()

    # Write the data to a png image called world.png the current directory
    mapnik.render_to_file(m,'sc.png', 'png')


def map_proj():
    import pyproj
    sc_proj='+proj=lcc +lat_1=32.5 +lat_2=34.83333333333334 +lat_0=31.83333333333333 +lon_0=-81 +x_0=609600 +y_0=0 +datum=NAD83 +units=ft +no_defs '
    projection=pyproj.Proj(sc_proj)
    print(projection(34.0, 81.0))

if __name__ == "__main__":
    map_proj()
