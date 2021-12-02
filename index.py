import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

from app import app
from apps import gvd, spm, nlse, dash_plot, modes_dsply, results


app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    html.Div(id='page-content')
])


@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
def display_page(pathname):
    if pathname == '/apps/gvd':
        return gvd.layout
    elif pathname == '/apps/spm':
         return spm.layout
    elif pathname == '/apps/nlse':
        return nlse.layout
    elif pathname == '/apps/dash_plot':
        return dash_plot.layout
    elif pathname == '/apps/modes_dsply':
        return modes_dsply.layout
    elif pathname == '/apps/results':
        return results.layout
    else:
        return '404'

if __name__ == '__main__':
    #app.run_server(host= '0.0.0.0',debug=False)#now you can run in other device using the same network. 
    app.run_server(host= '0.0.0.0',debug=True)#now you can run in other device using the same network. 
    #Changing the localhost (127.0.0.1:8050) to the IPV4 (abhängig von LAN- bzw. WLAN-Verbindung)
    #app.run_server(debug=True)
