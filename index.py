import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

from app import app
from apps import gvd, spm, nlse


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
    # elif pathname == '/apps/results':
    #     return results.layout
    else:
        return '404'

if __name__ == '__main__':
    app.run_server(debug=True)