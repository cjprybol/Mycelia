# # neo_import_dir = "/Users/cameronprybol/Library/Application Support/Neo4j Desktop/Application/relate-data/dbmss/dbms-8ab8baac-5dea-4137-bb24-e0b426447940/import"

# # uploading over API is slow for remote and local connections
# # Progress:   0%|▏                                        |  ETA: 8:13:39
# # upload_nodes_over_api(graph, ADDRESS=local_neo4j_bolt_address, PASSWORD=local_neo4j_password)
# # Progress:   0%|         

# # # push to Neo4J Aura
# # # run(`sudo touch /etc/neo4j/neo4j.conf`)
# # run(`sudo neo4j stop`)
# # # remote database needs to be running
# # # needs to be big enough
# # # leave off port from address
# # # run(`neo4j-admin push-to-cloud --overwrite --verbose --bolt-uri=$(ADDRESS) --username=$(USERNAME) --password=$(PASSWORD)`)
# # # run(`sudo neo4j-admin push-to-cloud --overwrite --verbose --dump-to "$(DIR)/test.db.dump" --bolt-uri=$(a) --username=$(USERNAME) --password=$(PASSWORD)`)
# # run(`sudo neo4j-admin push-to-cloud --overwrite --verbose --bolt-uri=$(a) --username=$(USERNAME) --password=$(PASSWORD)`)

# # cmd = "CREATE CONSTRAINT ON (t:Taxonomy) ASSERT t.tax_id IS UNIQUE"
# # @time run(cypher(address = local_address, username = USERNAME, password = local_password, database = DATABASE, cmd = cmd))

# # cmd = "CREATE CONSTRAINT ON (t:Taxonomy) ASSERT t.`scientific name` IS UNIQUE"
# # @time run(cypher(address = local_address, username = USERNAME, password = local_password, database = DATABASE, cmd = cmd))

# # cmd = "CREATE CONSTRAINT ON (t:Taxonomy) ASSERT t.identifier IS UNIQUE"
# # @time run(cypher(address = local_address, username = USERNAME, password = local_password, database = DATABASE, cmd = cmd))

# # parameters = ["$(n): row.$(n)" for n in filter(x -> x != "TYPE", names(node_table))]
# # parameters = "{" * join(parameters, ", ") * "}"

# # window_size = 10000
# # V = DataFrames.nrow(node_table)
# # windows = [i:min(i+window_size-1,V) for i in 1:window_size:V]
# # ProgressMeter.@showprogress for (i, w) in enumerate(windows)
# #     df_sub = node_table[w, :]
# #     f = "node$i.tsv"
# #     local_f_path = "$(temp_upload_dir)/$(f)"
# #     uCSV.write(local_f_path, df_sub, delim='\t')
# #     run(`chmod 777 $(local_f_path)`)
# #     f_url = "file:///$(local_f_path)"
# #     cmd =
# #     """
# #     LOAD CSV WITH HEADERS FROM '$(f_url)' AS row FIELDTERMINATOR '\t'
# #     CREATE (node:$(NODE_TYPE) $(parameters))
# #     """
# #     cmd = rstrip(replace(cmd, '\n' => ' '))
# #     cypher_cmd = Mycelia.cypher(address = local_address, username = USERNAME, password = local_password, database = DATABASE, cmd = cmd)
# #     run(cypher_cmd) 
# # end

# # ProgressMeter.@showprogress for (i, w) in enumerate(windows)
# #     df_sub = edge_table[w, :]
# #     f = "edge$i.tsv"
# #     local_f_path = "$(temp_upload_dir)/$(f)"
# #     uCSV.write(local_f_path, df_sub, delim='\t')
# #     run(`chmod 777 $(local_f_path)`)
# #     f_url = "file:///$(local_f_path)"
# #     cmd = 
# #     """
# #     LOAD CSV WITH HEADERS FROM '$(f_url)' AS row FIELDTERMINATOR '\t'
# #     MATCH (src:$(src_type) {identifier: row.src})
# #     MATCH (dst:$(dst_type) {identifier: row.dst})
# #     MERGE (src)-[p:$(edge_type)]->(dst)
# #     """
# #     cmd = rstrip(replace(cmd, '\n' => ' '))
# #     cypher_cmd = Mycelia.cypher(address = local_address, username = USERNAME, password = local_password, database = DATABASE, cmd = cmd)
# #     run(cypher_cmd) 
# # end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Upload all nodes from a MetaGraph to a Neo4j database, processing each unique node type separately.

# # Arguments
# - `graph`: A MetaGraph containing nodes to be uploaded
# - `address`: Neo4j server address (e.g., "bolt://localhost:7687")
# - `username`: Neo4j authentication username (default: "neo4j")
# - `password`: Neo4j authentication password
# - `format`: Data format for upload (default: "auto")
# - `database`: Target Neo4j database name (default: "neo4j")
# - `neo4j_import_directory`: Path to Neo4j's import directory for bulk loading
# """
# function upload_nodes_to_neo4j(;graph, address, username="neo4j", password, format="auto", database="neo4j", neo4j_import_directory)

#     node_types = unique(MetaGraphs.props(graph, v)[:TYPE] for v in Graphs.vertices(graph))
#     # node_type_strings = Mycelia.type_to_string.(node_types)

#     for node_type in node_types
#         @info "uploading node_type => $(Mycelia.type_to_string(node_type))..."
#         node_type_table = node_type_to_dataframe(node_type=node_type, graph=graph)
#         try
#             upload_node_table(table=node_type_table, address=address, password=password, neo4j_import_dir=neo4j_import_directory)
#         catch e
#             showerror(stdout, e)
#         end
#     end

#     @info "done!"
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Convert all nodes of a specific type in a MetaGraph to a DataFrame representation.

# # Arguments
# - `node_type`: The type of nodes to extract from the graph
# - `graph`: A MetaGraph containing the nodes

# # Returns
# A DataFrame where:
# - Each row represents a node of the specified type
# - Columns correspond to all unique properties found across nodes
# - Values are JSON-serialized strings for consistency

# # Notes
# - All values are normalized through JSON serialization
# - Dictionary values receive double JSON encoding
# - The TYPE column is converted using `type_to_string`
# """
# function node_type_to_dataframe(;node_type, graph)
#     node_type_indices = filter(v -> MetaGraphs.props(graph, v)[:TYPE] == node_type, Graphs.vertices(graph))
#     node_type_parameters = unique(reduce(vcat, map(v -> collect(keys(MetaGraphs.props(graph, v))), node_type_indices)))
#     node_type_table = DataFrames.DataFrame(Dict(p => [] for p in node_type_parameters))
#     for node_index in node_type_indices      
#         push!(node_type_table, MetaGraphs.props(graph, node_index))
#     end
#     # normalize
#     node_type_table[!, "TYPE"] = Mycelia.type_to_string.(node_type_table[!, "TYPE"])
#     for column in names(node_type_table)
#         T = Union{unique(typeof.(node_type_table[!, column]))...}
#         if T <: AbstractDict
#             node_type_table[!, column] = map(d -> JSON.json(string(JSON.json(d))), node_type_table[!, column])
#         else
#             node_type_table[!, column] = JSON.json.(string.(node_type_table[!, column]))
#         end
#     end  
#     return node_type_table
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Upload a DataFrame to Neo4j as nodes in batched windows.

# # Arguments
# - `table::DataFrame`: Input DataFrame where each row becomes a node. Must contain a "TYPE" column.
# - `address::String`: Neo4j server address (e.g. "bolt://localhost:7687")
# - `password::String`: Neo4j database password
# - `neo4j_import_dir::String`: Directory path accessible to Neo4j for importing files
# - `window_size::Int=1000`: Number of rows to process in each batch
# - `username::String="neo4j"`: Neo4j database username
# - `database::String="neo4j"`: Target Neo4j database name

# # Notes
# - All rows must have the same node type (specified in "TYPE" column)
# - Column names become node properties
# - Requires write permissions on neo4j_import_dir
# - Large tables are processed in batches of size window_size
# """
# function upload_node_table(;table, window_size=1000, address, password, username="neo4j", database="neo4j", neo4j_import_dir)
#     nrows = DataFrames.nrow(table)
#     windows = (i:min(i+window_size-1,nrows) for i in 1:window_size:nrows)

#     node_types = unique(table[!, "TYPE"])
#     @assert length(node_types) == 1
#     NODE_TYPE = Mycelia.type_to_string(first(node_types))
#     parameters = ["$(n): row.$(n)" for n in filter(x -> !(x in ["TYPE"]), names(table))]
#     parameters = "{" * join(parameters, ", ") * "}"

#     ProgressMeter.@showprogress for (i, window) in enumerate(windows)
#         df_sub = table[window, :]
#         f = "node$i.tsv"
#         local_f_path = "$(neo4j_import_dir)/$(f)"
#         uCSV.write(local_f_path, df_sub, delim='\t')
#         run(`chmod 777 $(local_f_path)`)
#         f_url = "file:///$(f)"
#         cmd =
#         """
#         LOAD CSV WITH HEADERS FROM '$(f_url)' AS row FIELDTERMINATOR '\t'
#         CREATE (:`$(NODE_TYPE)` $(parameters))
#         """
#         cmd = rstrip(replace(cmd, '\n' => ' '))
#         cypher_cmd = Mycelia.cypher(cmd, address = address, username = username, password = password, database = database)
#         run(cypher_cmd) 
#     end
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Upload a single node from a MetaGraph to a Neo4j database using the HTTP API.

# # Arguments
# - `graph`: MetaGraph containing the node to be uploaded
# - `v`: Vertex identifier in the graph
# - `ADDRESS`: Neo4j server address (e.g. "http://localhost:7474")
# - `USERNAME`: Neo4j authentication username (default: "neo4j")
# - `PASSWORD`: Neo4j authentication password
# - `DATABASE`: Target Neo4j database name (default: "neo4j")

# # Details
# Generates and executes a Cypher MERGE command using the node's properties. The node's :TYPE 
# and :identifier properties are used for node labeling, while other non-empty properties 
# are added as node properties.
# """
# function upload_node_over_api(graph, v; ADDRESS, USERNAME="neo4j", PASSWORD, DATABASE="neo4j")
#     node_type = MetaGraphs.props(graph, v)[:TYPE]
#     node_identifier = MetaGraphs.props(graph, v)[:identifier]
#     node_parameters = filter(x -> 
#             !(x[1] in (:TYPE, :identifier)) && 
#             !(ismissing(x[2]) || isempty(x[2])), 
#         MetaGraphs.props(graph, v))
#     params_string = join(["$(string(key)): \"$(string(value))\"" for (key, value) in node_parameters], ", ")
#     node_type_string = Mycelia.type_to_string(node_type)
#     node_identifier_string = string(node_identifier)
#     cmd = 
#     """
#     MERGE (`$(node_identifier_string)`:`$(node_type_string)` {$(params_string)})
#     """
#     cmd = strip(cmd)
#     cypher_cmd = Mycelia.cypher(cmd, address = ADDRESS, username = USERNAME, password = PASSWORD, database = DATABASE)
#     run(cypher_cmd)
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Uploads all nodes from the given graph to a specified API endpoint.

# # Arguments
# - `graph`: The graph containing the nodes to be uploaded.
# - `ADDRESS`: The API endpoint address.
# - `USERNAME`: The username for authentication (default: "neo4j").
# - `PASSWORD`: The password for authentication.
# - `DATABASE`: The database name (default: "neo4j").
# """
# function upload_nodes_over_api(graph; ADDRESS, USERNAME="neo4j", PASSWORD, DATABASE="neo4j")
#     ProgressMeter.@showprogress for v in Graphs.vertices(graph)
#         upload_node_over_api(graph, v, ADDRESS=ADDRESS, USERNAME=USERNAME, PASSWORD=PASSWORD, DATABASE=DATABASE)
#     end    
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Creates unique identifier constraints for each node type in a Neo4j database.

# # Arguments
# - `graph`: A MetaGraph containing nodes with TYPE properties
# - `address`: Neo4j server address
# - `username`: Neo4j username (default: "neo4j")
# - `password`: Neo4j password
# - `database`: Neo4j database name (default: "neo4j")

# # Details
# Extracts unique node types from the graph and creates Neo4j constraints ensuring
# each node of a given type has a unique identifier property.

# Failed constraint creation attempts are silently skipped.
# """
# function create_node_constraints(graph; address, username="neo4j", password, database="neo4j")
#     node_types = unique(MetaGraphs.props(graph, v)[:TYPE] for v in Graphs.vertices(graph))
#     node_type_strings = map(t -> Mycelia.type_to_string(t), node_types)
#     for t in node_type_strings
#         cmd = "CREATE CONSTRAINT ON (t:`$(t)`) ASSERT t.identifier IS UNIQUE"
#         try
#             cypher = Mycelia.cypher(cmd, address = address, username = username, password = password, database = database)
#             @show cypher
#             @time run(cypher)
#         catch
#             continue
#         end
#     end
# end

# # function batch_upload_edge_type_over_url_from_graph(src_type, dst_type, edge_type, graph, ADDRESS, USERNAME, PASSWORD, DATABASE; window_size=100)    
# #     src_nodes = filter(v -> MetaGraphs.get_prop(graph, v, :TYPE) == src_type, Graphs.vertices(graph))
# #     dst_nodes = filter(v -> MetaGraphs.get_prop(graph, v, :TYPE) == dst_type, Graphs.vertices(graph))
# #     edges_to_upload = []
# #     for src_node in src_nodes
# #         outneighbors = Graphs.outneighbors(graph, src_node)
# #         outneighbors = filter(outneighbor -> outneighbor in dst_nodes, outneighbors)
# #         for outneighbor in outneighbors
# #             this_edge = Graphs.Edge(src_node, outneighbor)
# #             @assert MetaGraphs.get_prop(graph, this_edge, :TYPE) == edge_type
# #             push!(edges_to_upload, this_edge)
# #         end
# #     end

# #     N = length(edges_to_upload)
# #     windows = [i:min(i+window_size-1,N) for i in 1:window_size:N]

# #     ProgressMeter.@showprogress for window in windows
# #         cmds = []
# #         for (i, e) in enumerate(edges_to_upload[window])
# #             edge_params = filter(p -> p != :TYPE, keys(MetaGraphs.props(graph, e)))
# #             params = ["$(string(param)):'$(MetaGraphs.get_prop(graph, e, param))'" for param in edge_params]
# #             joint_params = join(params, ", ")
# #             node_cmds = 
# #             """
# #             MERGE (src$(i):$(MetaGraphs.props(graph, e.src)[:TYPE]) {identifier: '$(MetaGraphs.props(graph, e.src)[:identifier])'})
# #             MERGE (dst$(i):$(MetaGraphs.props(graph, e.dst)[:TYPE]) {identifier: '$(MetaGraphs.props(graph, e.dst)[:identifier])'})
# #             """
# #             if !isempty(joint_params)
# #                 relationship_cmd = "MERGE (src$(i))-[r$(i):$(MetaGraphs.props(graph, e)[:TYPE]) {$(joint_params)}]->(dst$(i))"
# #             else
# #                 relationship_cmd = "MERGE (src$(i))-[r$(i):$(MetaGraphs.props(graph, e)[:TYPE])]->(dst$(i))"
# #             end
# #             cmd = node_cmds * relationship_cmd
# #             cmd = replace(cmd, '\n' => ' ')
# #             push!(cmds, cmd)
# #         end
# #         cmd = join(cmds, ' ')
# #         cypher_cmd = Mycelia.cypher(address = ADDRESS, username = USERNAME, password = PASSWORD, database = DATABASE, cmd = cmd)
# #         run(cypher_cmd)
# #     end    
# # end

# # function batch_upload_node_type_over_url_from_graph(node_type, graph, ADDRESS, USERNAME, PASSWORD, DATABASE, window_size=100)
# #     node_type_params = Set{Symbol}()
# #     vertices_of_type = [v for v in Graphs.vertices(graph) if (graph.vprops[v][:TYPE] == node_type)]

# #     V = length(vertices_of_type)
# #     windows = [i:min(i+window_size-1,V) for i in 1:window_size:V]

# #     ProgressMeter.@showprogress for window in windows
# #         cmds = []
# #         for (i, v) in enumerate(vertices_of_type[window])
# #             node_params = filter(p -> p != :TYPE, keys(MetaGraphs.props(graph, v)))
# #             params = ["$(string(param)):'$(MetaGraphs.get_prop(graph, v, param))'" for param in node_params]
# #             # params = ["$(string(param)):'$(escape_string(MetaGraphs.get_prop(graph, v, param)))'" for param in node_params]
# #             joint_params = join(params, ", ")
# #             cmd = "MERGE (node$(i):$(node_type) {$(joint_params)})"
# #             push!(cmds, cmd)
# #         end
# #         cmd = join(cmds, ' ')
# #         cypher_cmd = Mycelia.cypher(address = ADDRESS, username = USERNAME, password = PASSWORD, database = DATABASE, cmd = cmd)
# #         run(cypher_cmd)
# #     end    
# # end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Upload nodes of a specific type from a graph to a Neo4j database using MERGE operations.

# # Arguments
# - `node_type`: The type label for the nodes to upload
# - `graph`: Source MetaGraph containing the nodes
# - `ADDRESS`: Neo4j server address (e.g. "bolt://localhost:7687")
# - `PASSWORD`: Neo4j database password
# - `USERNAME="neo4j"`: Neo4j username (defaults to "neo4j")
# - `DATABASE="neo4j"`: Target Neo4j database name (defaults to "neo4j")
# - `window_size=100`: Number of nodes to upload in each batch (defaults to 100)

# # Details
# Performs batched uploads of nodes using Neo4j MERGE operations. Node properties are 
# automatically extracted from the graph vertex properties, excluding the 'TYPE' property.
# """
# function upload_node_type_over_url_from_graph(;node_type, graph, ADDRESS, USERNAME="neo4j", PASSWORD, DATABASE="neo4j", window_size=100)
#     node_type_params = Set{Symbol}()
#     vertices_of_type = [v for v in Graphs.vertices(graph) if (graph.vprops[v][:TYPE] == node_type)]

#     V = length(vertices_of_type)
#     windows = [i:min(i+window_size-1,V) for i in 1:window_size:V]

#     ProgressMeter.@showprogress for window in windows
#         cmds = []
#         for (i, v) in enumerate(vertices_of_type[window])
#             node_params = filter(p -> p != :TYPE, keys(MetaGraphs.props(graph, v)))
#             params = ["$(string(param)):'$(MetaGraphs.get_prop(graph, v, param))'" for param in node_params]
#             # params = ["$(string(param)):'$(escape_string(MetaGraphs.get_prop(graph, v, param)))'" for param in node_params]
#             joint_params = join(params, ", ")
#             cmd = "MERGE (node$(i):$(node_type) {$(joint_params)})"
#             push!(cmds, cmd)
#         end
#         cmd = join(cmds, ' ')
#         cypher_cmd = Mycelia.cypher(address = ADDRESS, username = USERNAME, password = PASSWORD, database = DATABASE, cmd = cmd)
#         run(cypher_cmd)
#     end    
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Upload edges of a specific type from a MetaGraph to a Neo4j database, batching uploads in windows.

# # Arguments
# - `src_type`: Type of source nodes to filter
# - `dst_type`: Type of destination nodes to filter  
# - `edge_type`: Type of edges to upload
# - `graph`: MetaGraph containing the nodes and edges
# - `ADDRESS`: Neo4j server URL
# - `USERNAME`: Neo4j username (default: "neo4j")
# - `PASSWORD`: Neo4j password
# - `DATABASE`: Neo4j database name (default: "neo4j")
# - `window_size`: Number of edges to upload in each batch (default: 100)

# # Details
# - Filters edges based on source, destination and edge types
# - Preserves all edge properties except :TYPE when uploading
# - Uses MERGE operations to avoid duplicate nodes/relationships
# - Uploads are performed in batches for better performance
# - Progress is shown via ProgressMeter

# # Returns
# Nothing
# """
# function upload_edge_type_over_url_from_graph(;src_type, dst_type, edge_type, graph, ADDRESS, USERNAME="neo4j", PASSWORD, DATABASE="neo4j", window_size=100)    
#     src_nodes = filter(v -> MetaGraphs.get_prop(graph, v, :TYPE) == src_type, Graphs.vertices(graph))
#     dst_nodes = filter(v -> MetaGraphs.get_prop(graph, v, :TYPE) == dst_type, Graphs.vertices(graph))
#     edges_to_upload = []
#     for src_node in src_nodes
#         outneighbors = Graphs.outneighbors(graph, src_node)
#         outneighbors = filter(outneighbor -> outneighbor in dst_nodes, outneighbors)
#         for outneighbor in outneighbors
#             this_edge = Graphs.Edge(src_node, outneighbor)
#             @assert MetaGraphs.get_prop(graph, this_edge, :TYPE) == edge_type
#             push!(edges_to_upload, this_edge)
#         end
#     end

#     N = length(edges_to_upload)
#     windows = (i:min(i+window_size-1,N) for i in 1:window_size:N)

#     ProgressMeter.@showprogress for window in windows
#         cmds = []
#         for (i, e) in enumerate(edges_to_upload[window])
#             edge_params = filter(p -> p != :TYPE, keys(MetaGraphs.props(graph, e)))
#             params = ["$(string(param)):'$(MetaGraphs.get_prop(graph, e, param))'" for param in edge_params]
#             joint_params = join(params, ", ")
#             node_cmds = 
#             """
#             MERGE (src$(i):$(MetaGraphs.props(graph, e.src)[:TYPE]) {identifier: '$(MetaGraphs.props(graph, e.src)[:identifier])'})
#             MERGE (dst$(i):$(MetaGraphs.props(graph, e.dst)[:TYPE]) {identifier: '$(MetaGraphs.props(graph, e.dst)[:identifier])'})
#             """
#             if !isempty(joint_params)
#                 relationship_cmd = "MERGE (src$(i))-[r$(i):$(MetaGraphs.props(graph, e)[:TYPE]) {$(joint_params)}]->(dst$(i))"
#             else
#                 relationship_cmd = "MERGE (src$(i))-[r$(i):$(MetaGraphs.props(graph, e)[:TYPE])]->(dst$(i))"
#             end
#             cmd = node_cmds * relationship_cmd
#             cmd = replace(cmd, '\n' => ' ')
#             push!(cmds, cmd)
#         end
#         cmd = join(cmds, ' ')
#         cypher_cmd = Mycelia.cypher(address = ADDRESS, username = USERNAME, password = PASSWORD, database = DATABASE, cmd = cmd)
#         run(cypher_cmd)
#     end    
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Constructs a command to execute Neo4j Cypher queries via cypher-shell.

# # Arguments
# - `cmd`: The Cypher query command to execute
# - `address::String="neo4j://localhost:7687"`: Neo4j server address
# - `username::String="neo4j"`: Neo4j authentication username
# - `password::String="password"`: Neo4j authentication password 
# - `format::String="auto"`: Output format (auto, verbose, or plain)
# - `database::String="neo4j"`: Target Neo4j database name

# # Returns
# - `Cmd`: A command object ready for execution
# """
# function cypher(cmd;
#     address="neo4j://localhost:7687",
#     username="neo4j",
#     password="password",
#     format="auto",
#     database="neo4j"
#     )
#     # cmd = `cypher-shell --address $(address) --username $(username) --password $(password) --format $(format) --database $(database) $(cmd)`
#     cmd = `/home/jovyan/.local/cypher-shell/cypher-shell --address $(address) --username $(username) --password $(password) --format $(format) --database $(database) $(cmd)`
# #    cmd = `/home/jupyter-cjprybol/software/cypher-shell/cypher-shell --address $(address) --username $(username) --password $(password) --format $(format) --database $(database) $(cmd)`
#     return cmd
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Lists all available Neo4j databases on the specified server.

# # Arguments
# - `address::String`: Neo4j server address (e.g. "neo4j://localhost:7687")
# - `username::String="neo4j"`: Neo4j authentication username
# - `password::String`: Neo4j authentication password

# # Returns
# - `DataFrame`: Contains database information with columns typically including:
#   - name: Database name
#   - address: Database address
#   - role: Database role (e.g., primary, secondary)
#   - status: Current status (e.g., online, offline)
#   - default: Boolean indicating if it's the default database
# """
# function list_databases(;address, username="neo4j", password)
#     cmd = "show databases"
#     database = "system"
#     cmd = cypher(cmd, address=address, username=username, password=password, database=database)
#     return DataFrames.DataFrame(uCSV.read(open(cmd), header=1, quotes='"', encodings=Dict("FALSE" => false, "TRUE" => true))...)
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Creates a new Neo4j database instance if it doesn't already exist.

# # Arguments
# - `database::String`: Name of the database to create
# - `address::String`: Neo4j server address (e.g. "neo4j://localhost:7687")
# - `username::String`: Neo4j authentication username (defaults to "neo4j")
# - `password::String`: Neo4j authentication password

# # Notes
# - Requires system database privileges to execute
# - Silently returns if database already exists
# - Temporarily switches to system database to perform creation
# """
# function create_database(;database, address, username="neo4j", password)
#     current_databases = list_databases(;address, username, password)
#     if database in current_databases[!, "name"]
#         return
#     else
#         f = run
#         cmd = "create database $(database)"
#         # switch database to system, so that we can create the user-specific database in the system
#         database = "system"
#         run(cypher(;address, username, password, database, cmd, f))
#     end
# end

# # """
# # $(DocStringExtensions.TYPEDSIGNATURES)

# # A short description of the function

# # ```jldoctest
# # julia> 1 + 1
# # 2
# # ```
# # """
# # function cypher(;address, username, password, database, cmd)
# #     return `cypher-shell --address $address --username $username --password $password --database $(database) --format auto $(cmd)`
# # end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Query Neo4j database to find all descendant taxonomic IDs for a given taxonomic ID.

# # Arguments
# - `tax_id`: Source taxonomic ID to find children for
# - `DATABASE_ID`: Neo4j database identifier (required)
# - `USERNAME`: Neo4j database username (default: "neo4j")
# - `PASSWORD`: Neo4j database password (required)

# # Returns
# `Vector{Int}`: Sorted array of unique child taxonomic IDs
# """
# function taxonomic_id_to_children(tax_id; DATABASE_ID, USERNAME="neo4j", PASSWORD)
#     DATABASE = "neo4j"
#     ADDRESS="neo4j+s://$(database_id).databases.neo4j.io:7687"

#     # NOTE! *, or 0 distance (e.g. [*0..2]) step range will include source node!!!!
#     cmd = "MATCH (n)<-[*]-(n2) WHERE n.tax_id IS NOT NULL AND n.tax_id = \"$(tax_id)\" RETURN DISTINCT n2.tax_id AS tax_id"
#     println(cmd)

#     cypher = cypher(cmd, address=ADDRESS, username = USERNAME, password = PASSWORD, database = DATABASE)
#     tax_ids = readlines(open(cypher))[2:end]
#     tax_ids = strip.(tax_ids, '"')
#     tax_ids = parse.(Int, tax_ids)
#     return unique(tax_ids)
# end
