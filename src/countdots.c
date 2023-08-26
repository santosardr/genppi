#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <ctype.h>
#include <stdbool.h>
#include <glib.h>

typedef struct
{
    char *node1;
    char *node2;
} Edge;

typedef struct
{
    GHashTable *node_set;
    GHashTable *edge_set;
    int num_nodes;
    long num_edges;
} Graph;

void add_node(Graph *graph, const char *node)
{
    if (!g_hash_table_contains(graph->node_set, node))
    {
        g_hash_table_add(graph->node_set, g_strdup(node));
        graph->num_nodes++;
    }
}

void add_edge(Graph *graph, const char *node1, const char *node2)
{
    
    // Sort the nodes before adding the edge
    char sorted_node1[256];
    char sorted_node2[256];
    if (strcmp(node1, node2) < 0)
    {
        strcpy(sorted_node1, node1);
        strcpy(sorted_node2, node2);
    }
    else
    {
        strcpy(sorted_node1, node2);
        strcpy(sorted_node2, node1);
    }
    // Create the edge key by concatenating the two node names
    char edge_key[512];
    snprintf(edge_key, sizeof(edge_key), "%s-%s", sorted_node1, sorted_node2);

    // Check if the edge already exists in the set
    if (g_hash_table_contains(graph->edge_set, edge_key)) {
        return;
    }
    // Add the edge to the set
    g_hash_table_add(graph->edge_set, g_strdup(edge_key));
    graph->num_edges++;
}

char *extractNodeName(const char *input)
{
    char *result = NULL;
    if (input)
    {
        // Remove any leading or trailing spaces
        const char *start = input;
        while (isspace(*start))
        {
            start++;
        }
        const char *end = input + strlen(input) - 1;
        while (end > start && isspace((unsigned char)*end))
        {
            end--;
        }

        // Check if the node name is enclosed in double quotes
        if (start[0] == '\"' && end[0] == '\"')
        {
            // Extract the node name without the double quotes
            result = (char *)malloc((end - start) * sizeof(char));
            strncpy(result, start + 1, end - start - 1);
            result[end - start - 1] = '\0';
        }
        else
        {
            // No double quotes, use the node name as is
            result = strdup(start);
        }
    }
    return result;
}

void count_nodes_edges(const char *dot_file, Graph *graph)
{
    FILE *file = fopen(dot_file, "r");
    if (file == NULL)
    {
        printf("Failed to open file: %s\n", dot_file);
        exit(1);
    }

    // Create the  sets
    graph->edge_set = g_hash_table_new(g_str_hash, g_str_equal);
    graph->node_set = g_hash_table_new(g_str_hash, g_str_equal);

    char line[256];
    while (fgets(line, sizeof(line), file))
    {
        char *node1 = NULL;
        char *node2 = NULL;

        // Check if the line contains an edge connection
        if (strstr(line, " -- "))
        {
            // Split the line by " -- "
            char *token = strtok(line, " -- ");
            node1 = extractNodeName(token);
            token = strtok(NULL, " -- ");
            node2 = extractNodeName(token);
        }

        if (node1 != NULL && node2 != NULL)
        {
            add_edge(graph, node1, node2);
            add_node(graph, node1);
            add_node(graph, node2);
        }
    }

    fclose(file);
}

void process_file(const char *dot_file, int print_labels)
{
    Graph *graph = (Graph *)malloc(sizeof(Graph));
    if (graph == NULL)
    {
        printf("Memory allocation failed\n");
        return;
    }
    graph->num_nodes = 0;
    graph->num_edges = 0;
    graph->edge_set = NULL;

    count_nodes_edges(dot_file, graph);

    if (print_labels)
    {
        printf("File\tNodes\tEdges\n");
    }

    char *filename = strrchr(dot_file, '/');
    if (filename == NULL)
    {
        filename = (char *)dot_file;
    }
    else
    {
        filename++;
    }

    printf("%s\t%d\t%ld\n", filename, graph->num_nodes, graph->num_edges);
    free(graph);
}

void process_folder(const char *folder)
{
    int print_labels = 1;
    int total_nodes = 0;
    int total_edges = 0;

    char command[256];
    sprintf(command, "find %s -name '*.dot'", folder);
    FILE *pipe = popen(command, "r");
    if (pipe == NULL)
    {
        printf("Failed to execute command: %s\n", command);
        exit(1);
    }
    char dot_file[256];
    while (fgets(dot_file, sizeof(dot_file), pipe))
    {
        dot_file[strcspn(dot_file, "\n")] = 0;

        process_file(dot_file, print_labels);
        print_labels = 0;

        Graph *graph = (Graph *)malloc(sizeof(Graph));
        if (graph == NULL)
        {
            printf("Memory allocation failed\n");
            return;
        }

        graph->num_nodes = 0;
        graph->num_edges = 0;
        graph->edge_set = NULL;
        count_nodes_edges(dot_file, graph);

        total_nodes += graph->num_nodes;
        total_edges += graph->num_edges;

        free(graph);
    }

    pclose(pipe);

    if (total_nodes > 0 || total_edges > 0)
    {
        printf("\nSummary - Total nodes: %d, Total edges: %d\n", total_nodes, total_edges);
    }
}
int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        printf("Usage: ./countdots <dot_file_or_folder>\n");
        return 1;
    }

    char *path = argv[1];

    if (access(path, F_OK) != 0)
    {
        printf("Invalid file or folder path.\n");
        return 1;
    }

    if (access(path, R_OK) != 0)
    {
        printf("Permission denied: %s\n", path);
        return 1;
    }
    struct stat path_stat;
    if (stat(path, &path_stat) == 0)
    {
        if (S_ISREG(path_stat.st_mode))
        {
            process_file(path, 1);
        }
        else if (S_ISDIR(path_stat.st_mode))
        {
            process_folder(path);
        }
        else
        {
            printf("Invalid file or folder path.\n");
            return 1;
        }
    }
    else
    {
        printf("Failed to get file or directory status.\n");
        return 1;
    }
}
